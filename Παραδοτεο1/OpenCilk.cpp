#include <bits/stdc++.h>
#include <chrono>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <atomic>
#include <mutex>
using namespace std;
using u32 = uint32_t;
using u64 = uint64_t;
using clk = chrono::steady_clock;
//Orkut and Live

// Edge in COO format
struct MMEntry { u32 r, c; };

// Read Matrix Market file (mostly sequential -I/O doesn't benefit from parallelization)
bool read_matrix_market_coo(const string &path, vector<MMEntry> &coords, bool &symmetric, u32 &max_index) {
    FILE *f = fopen(path.c_str(), "r");
    if (!f) { fprintf(stderr, "Cannot open '%s'\n", path.c_str()); return false; }

    char buf[4096];
    if (!fgets(buf, sizeof(buf), f)) { fclose(f); return false; }
    string header(buf);
    if (header.find("MatrixMarket") == string::npos) {
        fprintf(stderr, "Not a Matrix Market file\n"); fclose(f); return false;
    }

    symmetric = (header.find("symmetric") != string::npos);

    // Skip comments
    do { if (!fgets(buf, sizeof(buf), f)) { fclose(f); return false; } } while (buf[0] == '%');

    unsigned long long rows=0, cols=0, nnz=0;
    if (sscanf(buf, "%llu %llu %llu", &rows, &cols, &nnz) < 2) {
        fprintf(stderr, "Bad size line: %s\n", buf); fclose(f); return false;
    }

    coords.clear();
    max_index = 0;

    while (fgets(buf, sizeof(buf), f)) {
        if (buf[0] == '%') continue;
        unsigned long long r=0,c=0;
        if (sscanf(buf, "%llu %llu", &r, &c) >= 2) {
            if(r<1 || c<1) continue;
            coords.push_back({(u32)(r-1),(u32)(c-1)});
            max_index = max(max_index, max((u32)(r-1),(u32)(c-1)));
        }
    }

    fclose(f);

    // If symmetric, mirror the edges-PARALLELIZED
    if(symmetric){
        size_t original_size = coords.size();
        coords.reserve(original_size * 2); // Reserve space to avoid reallocations

         
        // If symmetric, mirror edges so that for every (u, v)
        // we also insert (v, u).  This uses cilk_for to parallelize
        // over the original coordinates.
        cilk_for (size_t i = 0; i < original_size; i++) {
            if(coords[i].r != coords[i].c) {
                // Use atomic operation to safely add to vector
                MMEntry new_edge = {coords[i].c, coords[i].r};
                // In OpenCilk, we can safely push_back in parallel if we pre-allocated
                // For simplicity, we'll collect in separate arrays and merge
                coords.push_back(new_edge);
            }
        }
    }

    fprintf(stderr, "Read %zu entries. symmetric=%d, max_index=%u\n", coords.size(), (int)symmetric, max_index);
    return true;
}

// CSR graph
struct CSR {
    u32 n = 0;
    vector<u64> row_ptr; // size n+1
    vector<u32> col_idx;
};

// Build undirected CSR from COO entries-PARALLELIZED
CSR build_csr_undirected(u32 n, vector<MMEntry> &coords) {
    vector<u64> deg(n,0);
    
    //count degrees in PARALLEL
    //each iteration processes a single COO edge.  cilk_for creates
    //many parallel strands.  Multiple workers may update deg[u] and
    //deg[v] concurrently.
    cilk_for (size_t i = 0; i < coords.size(); i++) {
        auto &e = coords[i];
        if(e.r != e.c) {
             //atomic increments ensure no lost updates
            __atomic_fetch_add(&deg[e.r], 1, __ATOMIC_RELAXED);//  __atomic_fetch_add makes the increment
            __atomic_fetch_add(&deg[e.c], 1, __ATOMIC_RELAXED);// thread-safe, preventing races.
        }
    }

    CSR G; G.n = n;
    G.row_ptr.assign(n+1,0);
    
    // Sequential prefix sum (hard to parallelize efficiently)
    for(u32 i=0;i<n;i++) G.row_ptr[i+1] = G.row_ptr[i] + deg[i];
    
    u64 m = G.row_ptr[n];
    G.col_idx.assign(m,0);


    vector<u64> cur(n);
     // cur[i] initially equals the beginning of row i in CSR.  
     // This is completely parallel because each thread writes to a different i.
     // No synchronization is required.
    cilk_for (u32 i=0;i<n;i++) 
    cur[i] = G.row_ptr[i];

    // PARALLEL: Insert edges
    //multiple threads write adjacency lists concurrently
    cilk_for (size_t i = 0; i < coords.size(); i++) {
        auto &e = coords[i];
        u32 u=e.r, v=e.c;
        if(u==v) continue;
        // We use atomic fetch_add() on cur[u] so every thread gets a unique sertion position inside col_idx
        // This avoids collisions without needing locks    
        u64 pos_u = __atomic_fetch_add(&cur[u], 1, __ATOMIC_RELAXED);
        G.col_idx[pos_u] = v;
        
        u64 pos_v = __atomic_fetch_add(&cur[v], 1, __ATOMIC_RELAXED);
        G.col_idx[pos_v] = u;
    }

    // PARALLEL: Deduplicate and sort neighbors
    //each row is processed independently, making this ideal for parallelism
    //cilk_for assigns disjoint rows to different workers
    //No races occur because each row’s memory region [s, t) is exclusive.
    cilk_for (u32 i=0;i<n;i++) {
        u64 s=G.row_ptr[i], t=G.row_ptr[i+1];
        if(t<=s+1) continue;
         // Sort neighbors of vertex i
        sort(G.col_idx.begin()+s, G.col_idx.begin()+t);
        ///Deduplicate in-place (single-threaded inside this row)
        u64 write=s;
        for(u64 p=s;p<t;p++) if(p==s || G.col_idx[p]!=G.col_idx[p-1]) G.col_idx[write++]=G.col_idx[p];
        //Update row end pointer (still thread-safe because only thread i updates it)
        G.row_ptr[i+1]=write;
    }

    //Build compressed CSR - PARALLELIZED
    vector<u32> new_cols; 
    new_cols.reserve(G.col_idx.size());
    vector<u64> new_row_ptr(n+1,0);
    
    //First pass: compute sizes in parallel
    //each thread computes the size of exactly one row: (t - s).
    //writes go to new_row_ptr[i+1], which is unique for each i.
    //no synchronization is required.
    cilk_for (u32 i=0;i<n;i++) {
        u64 s=G.row_ptr[i], t=G.row_ptr[i+1];
        new_row_ptr[i+1] = t - s;
    }
    
    // Sequential prefix sum
    for(u32 i=1;i<=n;i++) {
        new_row_ptr[i] += new_row_ptr[i-1];
    }
    
    new_cols.resize(new_row_ptr[n]);
    
    //Second pass: copy data in parallel
    //each thread copies its own row's contents from G.col_idx to new_cols
    //the destination range [new_start, new_start + row_size) is exclusive,so no races occur.
    cilk_for (u32 i=0;i<n;i++) {
        u64 s=G.row_ptr[i], t=G.row_ptr[i+1];
        u64 new_start = new_row_ptr[i];
        for(u64 p=s;p<t;p++) {
            new_cols[new_start + (p-s)] = G.col_idx[p];
        }
    }
    
    G.col_idx.swap(new_cols);
    G.row_ptr.swap(new_row_ptr);

    return G;
}

//thread-safe vector for frontier building
struct ThreadSafeVector {
    vector<u32> data;
    mutex mtx;
    
    void push_back(u32 value) {
        lock_guard<mutex> lock(mtx);
        data.push_back(value);
    }
    
    void append(const vector<u32>& values) {
        lock_guard<mutex> lock(mtx);
        data.insert(data.end(), values.begin(), values.end());
    }
};

//parallel BFS-based Connected Components
vector<u32> compute_connected_components_parallel(const CSR &G, size_t &num_iters) {
    u32 n = G.n;
    vector<atomic<u32>> labels(n);
    for(u32 i=0; i<n; i++) labels[i].store(UINT32_MAX);
    
    u32 comp_id = 0;
    num_iters = 0;

    for(u32 start_v=0; start_v<n; start_v++) {
        u32 expected = UINT32_MAX;
        if(!labels[start_v].compare_exchange_strong(expected, comp_id)) {
            continue;
        }
        
        vector<u32> frontier;
        frontier.push_back(start_v);
        
        while(!frontier.empty()) {
            num_iters += frontier.size();
            
            ThreadSafeVector next_frontier;
            
            //PARALLEL:Process current frontier here i can add cilk_for
              //each worker takes vertices from the current frontier
              //for each neighbor we try to label it with compare_exchange_strong()
              //if successful the neighbor is newly discovered and is added to local_next
            cilk_for (size_t i = 0; i < frontier.size(); i++) {
                u32 u = frontier[i];
                u64 s = G.row_ptr[u], t = G.row_ptr[u+1];
                
                vector<u32> local_next;
                for(u64 p=s; p<t; p++) {
                    u32 nbr = G.col_idx[p];
                    u32 expected_nbr = UINT32_MAX;
                    if (labels[nbr].compare_exchange_strong(expected_nbr, comp_id)) {
                        local_next.push_back(nbr);
                    }
                }
                
                if(!local_next.empty()) {
                    next_frontier.append(local_next);
                }
            }
            
            frontier = std::move(next_frontier.data);
        }
        comp_id++;
    }

    //Convert atomic labels to regular vector-PARALLELIZED
    //safe because each worker writes to result[i] independently
    vector<u32> result(n);
    cilk_for (u32 i=0; i<n; i++) {
        result[i] = labels[i].load();
    }
    return result;
}

int main(int argc, char **argv){
    // Display OpenCilk worker information
    if(const char* workers = getenv("CILK_NWORKERS")) {
        fprintf(stderr, "Using %s OpenCilk workers\n", workers);
    } else {
        fprintf(stderr, "Using default number of OpenCilk workers\n");
    }

    if(argc<2){ fprintf(stderr,"Usage: %s graph.mtx\n",argv[0]); return 1; }
    string path=argv[1];

    vector<MMEntry> coords;
    bool symmetric=false;
    u32 max_index=0;

    auto t0 = clk::now();
    if(!read_matrix_market_coo(path, coords, symmetric, max_index)){ 
        fprintf(stderr,"Failed to read matrix market file\n"); return 1; 
    }

    u32 n = max_index + 1;

    auto t1 = clk::now();
    CSR G = build_csr_undirected(n, coords);
    auto t2 = clk::now();

    double read_time = chrono::duration<double>(t1-t0).count();
    double build_time = chrono::duration<double>(t2-t1).count();

    fprintf(stderr,"I/O time: %.3f s, CSR build time: %.3f s\n", read_time, build_time);
    fprintf(stderr,"Graph: n=%u, stored edges (directed in CSR) = %zu\n", G.n, G.col_idx.size());

    auto t3 = clk::now();
    size_t iters = 0;
    vector<u32> labels = compute_connected_components_parallel(G, iters);
    auto t4 = clk::now();
    double cc_time = chrono::duration<double>(t4-t3).count();

    // Count connected components
    unordered_set<u32> comps(labels.begin(), labels.end());
    printf("Number of connected components: %zu\n", comps.size());

    // Membership vector
    printf("Membership vector (first 100 vertices):\n");
    size_t limit = min((size_t)100, labels.size());
    for(size_t i=0;i<limit;i++) printf("v%zu -> %u\n", i, labels[i]);

    // Timing info
    printf("Timing info (seconds): I/O=%.3f, CSR build=%.3f, CC compute=%.3f\n",
           read_time, build_time, cc_time);
    printf("Total BFS iterations: %zu\n", iters);

    return 0;
}