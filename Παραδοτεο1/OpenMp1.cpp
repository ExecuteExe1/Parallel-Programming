#include <bits/stdc++.h>
#include <chrono>
#include <omp.h>
using namespace std;
using u32 = uint32_t;
using u64 = uint64_t;
using clk = chrono::steady_clock;
//For orkut n Live
// Edge in COO format
struct MMEntry { u32 r, c; };

// Read Matrix Market file (unchanged from sequential)
bool read_matrix_market_coo(const string &path, vector<MMEntry> &coords, bool &symmetric, u32 &max_index) {
    FILE *f = fopen(path.c_str(), "r");  //open file
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

    // If symmetric, mirror the edges same as the sequential
    if(symmetric){
        size_t original_size = coords.size();
        for(size_t i=0;i<original_size;i++){
            if(coords[i].r != coords[i].c)
                coords.push_back({coords[i].c, coords[i].r});
        }
    }

    fprintf(stderr, "Read %zu entries. symmetric=%d, max_index=%u\n", coords.size(), (int)symmetric, max_index);
    return true;
}

// CSR graph  same on the constructor
struct CSR {
    u32 n = 0;
    vector<u64> row_ptr; // size n+1
    vector<u32> col_idx;
};

//This is where the parallelism starts with OpenMp
CSR build_csr_undirected(u32 n, vector<MMEntry> &coords) {
    vector<u64> deg(n,0);
    
    //here OpenMP spawns a team of threads and each thread executes the code in the bracets { {}
    #pragma omp parallel //shared variables stay shared unless they somehow become private or thread-local
    {
        vector<u64> deg_private(n,0); //line runs once per thread
                                    //Basicaly each thread keeps its own counts to avoid conflicting writes to deg
        //this part basically distributes loop iterations across threads
        #pragma omp for nowait //this splits the loop interation among threads,how it works is thread[0] takes a chunk,thread[1] takes another etc etc
        for (size_t i = 0; i < coords.size(); i++) {
            auto &e = coords[i];
            if(e.r != e.c) {
                deg_private[e.r]++;
                deg_private[e.c]++;
            }
        } //we also avoid race conditions due to deg_private
        
        #pragma omp critical  //one thread at a time to prevent concurrent  writes on the shared deg[]
        {
            for(u32 i=0; i<n; i++) {  //each thread enters the critical section after finishing its loop work
                deg[i] += deg_private[i];  //then merges its private partial results into the shared array
            }
        }
    }

    CSR G; G.n = n;
    G.row_ptr.assign(n+1,0);
    for(u32 i=0;i<n;i++) G.row_ptr[i+1] = G.row_ptr[i] + deg[i];
    u64 m = G.row_ptr[n];
    G.col_idx.assign(m,0);

    vector<u64> cur(n);


    #pragma omp parallel for //creates a parallel loop
    for(u32 i=0;i<n;i++) cur[i] = G.row_ptr[i];  //each thread writes to disjoint elements of cur[i], so no synchronization is needed.

    // Parallel edge insertion,each thread gets a subset of edges to process.
    #pragma omp parallel for //the threads run independently except for controlled atomic operations inside.
    for(size_t idx=0; idx<coords.size(); idx++){
        auto &e = coords[idx];
        u32 u=e.r, v=e.c;
        if(u==v) continue;
        
        u64 pos_u, pos_v;
        #pragma omp atomic capture //ensures that the read (pos_u = cur[u]) and the increment (cur[u]++) happen atomically as one indivisible operation,which means that it prevents race conditions where two threads would try to insert into the same adjacency list.
        { pos_u = cur[u]; cur[u]++; } //without this atomic capture, two threads could overwrite each other.
        
        #pragma omp atomic capture  
        { pos_v = cur[v]; cur[v]++; }
        
        G.col_idx[pos_u] = v;
        G.col_idx[pos_v] = u;
    }

    // Parallel deduplication and sorting aka runs loop over all vertices in parallel,so each iteration handles the adjacency list of one vertex only--> no data races
    #pragma omp parallel for schedule(dynamic, 100) //each thread grabs chunks of 100 vertices at a time
    for(u32 i=0;i<n;i++){
        u64 s=G.row_ptr[i], t=G.row_ptr[i+1];
        if(t<=s+1) continue;
        sort(G.col_idx.begin()+s, G.col_idx.begin()+t);//Each thread grabs chunks of 100 vertices at a time
        u64 write=s;
        for(u64 p=s;p<t;p++) if(p==s || G.col_idx[p]!=G.col_idx[p-1]) G.col_idx[write++]=G.col_idx[p];
        G.row_ptr[i+1]=write;
    }

    // Sequential compaction (could be parallelized but might not be worth it)
    vector<u32> new_cols; new_cols.reserve(G.col_idx.size());
    vector<u64> new_row_ptr(n+1,0);
    new_row_ptr[0]=0;
    for(u32 i=0;i<n;i++){
        u64 s=G.row_ptr[i], t=G.row_ptr[i+1];
        for(u64 p=s;p<t;p++) new_cols.push_back(G.col_idx[p]);
        new_row_ptr[i+1] = new_cols.size();
    }
    G.col_idx.swap(new_cols);
    G.row_ptr.swap(new_row_ptr);

    return G;
}

// Parallel BFS-based Connected Components using label propagation
vector<u32> compute_connected_components(const CSR &G, size_t &num_iters){
    u32 n = G.n;
    vector<u32> labels(n);
    
    // Initialize labels with vertex IDs
    #pragma omp parallel for //assigns each vertex’s label in parallel.
    for(u32 i=0; i<n; i++) { //no races because each thread writes to different elements
        labels[i] = i;
    }

    bool changed = true;
    num_iters = 0;
    
    while(changed) {
        changed = false;
        num_iters++;
        
        #pragma omp parallel for schedule(dynamic, 100) //each thread processes a subset of vertices
        for(u32 u=0; u<n; u++) {
            u64 s = G.row_ptr[u], t = G.row_ptr[u+1];
            u32 current_label = labels[u];
            u32 min_label = current_label;
            
            // Find minimum label among neighbors
            for(u64 p=s; p<t; p++) {
                u32 nbr = G.col_idx[p];
                if(labels[nbr] < min_label) {
                    min_label = labels[nbr];
                }
            }
            
            // Update if we found a smaller label
            if(min_label < current_label) {
                labels[u] = min_label;
                #pragma omp atomic write
                changed = true;
            }
        }
        
        // Synchronize - additional iteration to propagate changes
        #pragma omp parallel for schedule(dynamic, 100)
        for(u32 u=0; u<n; u++) {
            u64 s = G.row_ptr[u], t = G.row_ptr[u+1];
            u32 current_label = labels[u];
            
            for(u64 p=s; p<t; p++) {
                u32 nbr = G.col_idx[p];
                if(labels[nbr] < current_label) {
                    current_label = labels[nbr];
                }
            }
            
            if(current_label < labels[u]) {
                labels[u] = current_label;
            }
        }
    }

    // Relabel to consecutive component IDs
    vector<u32> unique_labels(labels.begin(), labels.end());
    sort(unique_labels.begin(), unique_labels.end());
    unique_labels.erase(unique(unique_labels.begin(), unique_labels.end()), unique_labels.end());
    
    unordered_map<u32, u32> label_map;
    for(u32 i=0; i<unique_labels.size(); i++) {
        label_map[unique_labels[i]] = i;
    }
    
    #pragma omp parallel for //parallel rewrite of all labels,safe because each thread writes to a different index
    for(u32 i=0; i<n; i++) {
        labels[i] = label_map[labels[i]];
    }
    
    return labels;
}

int main(int argc, char **argv){
    if(argc<2){ fprintf(stderr,"Usage: %s graph.mtx\n",argv[0]); return 1; }
    string path=argv[1];

    // Set number of threads (optional)
    if(argc>2) {
        omp_set_num_threads(atoi(argv[2]));
    }

    vector<MMEntry> coords;
    bool symmetric=false;
    u32 max_index=0;
    //the rest is basically the same
    auto t0 = clk::now();
    if(!read_matrix_market_coo(path, coords, symmetric, max_index)){ fprintf(stderr,"Failed to read matrix market file\n"); return 1; }

    u32 n = max_index + 1;

    auto t1 = clk::now();
    CSR G = build_csr_undirected(n, coords);
    auto t2 = clk::now();

    double read_time = chrono::duration<double>(t1-t0).count();
    double build_time = chrono::duration<double>(t2-t1).count();

    fprintf(stderr,"I/O time: %.3f s, CSR build time: %.3f s\n", read_time, build_time);
    fprintf(stderr,"Graph: n=%u, stored edges (directed in CSR) = %zu\n", G.n, G.col_idx.size());
    fprintf(stderr,"Using %d OpenMP threads\n", omp_get_max_threads());

    auto t3 = clk::now();
    size_t iters = 0;
    vector<u32> labels = compute_connected_components(G, iters);
    auto t4 = clk::now();
    double cc_time = chrono::duration<double>(t4-t3).count();

    // Count connected components
    unordered_set<u32> comps(labels.begin(), labels.end());
    printf("Number of connected components: %zu\n", comps.size());
    printf("Number of iterations: %zu\n", iters);

    // Membership vector
    printf("Membership vector (label[v] for each vertex):\n");
    size_t limit = min((size_t)1000, labels.size());
    for(size_t i=0;i<limit;i++) printf("v%zu -> %u\n", i, labels[i]);

    // Timing info
    printf("Timing info (seconds): I/O=%.3f, CSR build=%.3f, CC compute=%.3f\n",
           read_time, build_time, cc_time);

    return 0;
}