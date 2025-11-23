#include <bits/stdc++.h>
#include <chrono>
#include <pthread.h>
using namespace std;
using u32 = uint32_t;
using u64 = uint64_t;
using clk = chrono::steady_clock;

//Sequencial Part
struct MMEntry { u32 r, c; };

//Read Matrix Market file
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

//CSR Structure
struct CSR {
    u32 n = 0;
    vector<u64> row_ptr;
    vector<u32> col_idx;
};

// Threaded CSR Structure
struct ThreadData {
    size_t thread_id;
    size_t num_threads;
    vector<MMEntry>* coords;
    vector<u64>* deg;
    vector<u64>* cur;
    CSR* G;
    u32 n;
};

//Degree counting
void* degree_count_thread(void* arg){//creating a custom thread structure,thread id,number of threads,pointer to the coords vector (edge list)pointer to deg array (global degree array)
    ThreadData* td = (ThreadData*)arg; //n=number of vertices
    size_t start = td->thread_id * td->coords->size() / td->num_threads; 
      //this splits the edge list among threads evenly


n = number of vertices
    size_t end   = (td->thread_id+1) * td->coords->size() / td->num_threads;

    vector<u64> local_deg(td->n,0);   //each thread has its own local degree array,avoids atomic operations inside the main loop (very slow) and each thread counts degree contributions locally

    for(size_t i=start;i<end;i++){
        auto &e = (*td->coords)[i];
        if(e.r != e.c){
            local_deg[e.r]++;
            local_deg[e.c]++;
        }
    }

    for(size_t i=0;i<td->n;i++)
    //this loop accumulates the thread’s local degrees into the global deg array
        __sync_fetch_and_add(&(*td->deg)[i], local_deg[i]); //__sync_fetch_and_add() is a GCC built-in atomic primitive
    return nullptr;
}

//Edge insertion 
void* edge_insert_thread(void* arg){//each thread inserts edges from its assigned chunk of the COO edge list into the CSR col_idx array in a thread-safe manner using atomic counters.
    ThreadData* td = (ThreadData*)arg;
    size_t start = td->thread_id * td->coords->size() / td->num_threads;
    size_t end   = (td->thread_id+1) * td->coords->size() / td->num_threads; //each thread gets a contiguous chunk of the edges

    for(size_t idx=start; idx<end; idx++){
        auto &e = (*td->coords)[idx];
        u32 u=e.r, v=e.c;
        if(u==v) continue;

        u64 pos_u = __sync_fetch_and_add(&(*td->cur)[u],1);
        u64 pos_v = __sync_fetch_and_add(&(*td->cur)[v],1);

        (*td->G).col_idx[pos_u] = v;
        (*td->G).col_idx[pos_v] = u;
    }
    return nullptr;
}

//Per-row sorting + deduplication
void* sort_dedup_thread(void* arg){
    ThreadData* td = (ThreadData*)arg;
    size_t start = td->thread_id * td->n / td->num_threads;
    size_t end   = (td->thread_id+1) * td->n / td->num_threads;

    for(u32 i=start;i<end;i++){
        u64 s = td->G->row_ptr[i], t = td->G->row_ptr[i+1];
        if(t<=s+1) continue;
        sort(td->G->col_idx.begin()+s, td->G->col_idx.begin()+t);
        u64 write = s;
        for(u64 p=s;p<t;p++)
            if(p==s || td->G->col_idx[p]!=td->G->col_idx[p-1]) td->G->col_idx[write++] = td->G->col_idx[p];
        td->G->row_ptr[i+1] = write;
    }
    return nullptr;
}

CSR build_csr_undirected_pthread(u32 n, vector<MMEntry> &coords, int num_threads){
    vector<u64> deg(n,0);
    pthread_t threads[num_threads];
    ThreadData td[num_threads];

    //Degree counting
    for(int i=0;i<num_threads;i++){
        td[i] = { (size_t)i, (size_t)num_threads, &coords, &deg, nullptr, nullptr, n };
        pthread_create(&threads[i], nullptr, degree_count_thread, &td[i]);
    }
    for(int i=0;i<num_threads;i++) pthread_join(threads[i], nullptr);

    //Build row_ptr
    CSR G; G.n = n;
    G.row_ptr.resize(n+1,0);
    for(u32 i=0;i<n;i++) G.row_ptr[i+1] = G.row_ptr[i] + deg[i];

    G.col_idx.resize(G.row_ptr[n],0);
    vector<u64> cur(n);
    for(u32 i=0;i<n;i++) cur[i] = G.row_ptr[i];

    //Edge insertion
    for(int i=0;i<num_threads;i++){
        td[i] = { (size_t)i, (size_t)num_threads, &coords, nullptr, &cur, &G, n };
        pthread_create(&threads[i], nullptr, edge_insert_thread, &td[i]);
    }
    for(int i=0;i<num_threads;i++) pthread_join(threads[i], nullptr);

    //Sorting + deduplication
    for(int i=0;i<num_threads;i++){
        td[i] = { (size_t)i, (size_t)num_threads, nullptr, nullptr, nullptr, &G, n };
        pthread_create(&threads[i], nullptr, sort_dedup_thread, &td[i]);
    }
    for(int i=0;i<num_threads;i++) pthread_join(threads[i], nullptr);

    //Compaction
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

//Sequential Connected Components
vector<u32> compute_connected_components(const CSR &G, size_t &num_iters){
    u32 n = G.n;
    vector<u32> labels(n);
    for(u32 i=0;i<n;i++) labels[i] = i;

    bool changed = true;
    num_iters = 0;
    while(changed){
        changed=false; num_iters++;
        for(u32 u=0; u<n; u++){
            u64 s=G.row_ptr[u], t=G.row_ptr[u+1];
            u32 min_label = labels[u];
            for(u64 p=s;p<t;p++) min_label = min(min_label, labels[G.col_idx[p]]);
            if(min_label < labels[u]){
                labels[u] = min_label;
                changed = true;
            }
        }
        for(u32 u=0; u<n; u++){
            u64 s=G.row_ptr[u], t=G.row_ptr[u+1];
            u32 cur = labels[u];
            for(u64 p=s;p<t;p++) cur = min(cur, labels[G.col_idx[p]]);
            labels[u] = cur;
        }
    }

    vector<u32> unique_labels(labels.begin(), labels.end());
    sort(unique_labels.begin(), unique_labels.end());
    unique_labels.erase(unique(unique_labels.begin(), unique_labels.end()), unique_labels.end());
    unordered_map<u32,u32> label_map;
    for(u32 i=0;i<unique_labels.size();i++) label_map[unique_labels[i]] = i;
    for(u32 i=0;i<n;i++) labels[i] = label_map[labels[i]];
    return labels;
}

/ Main 
int main(int argc, char **argv){
    if(argc<2){ fprintf(stderr,"Usage: %s graph.mtx [num_threads]\n",argv[0]); return 1; }
    string path = argv[1];
    int num_threads = 4;
    if(argc>2) num_threads = atoi(argv[2]);

    vector<MMEntry> coords;
    bool symmetric=false;
    u32 max_index=0;

    auto t0 = clk::now();
    if(!read_matrix_market_coo(path, coords, symmetric, max_index)){ fprintf(stderr,"Failed to read matrix\n"); return 1; }
    u32 n = max_index + 1;

    auto t1 = clk::now();
    CSR G = build_csr_undirected_pthread(n, coords, num_threads);
    auto t2 = clk::now();

    double read_time = chrono::duration<double>(t1-t0).count();
    double build_time = chrono::duration<double>(t2-t1).count();

    fprintf(stderr,"I/O time: %.3f s, CSR build time: %.3f s\n", read_time, build_time);
    fprintf(stderr,"Graph: n=%u, stored edges (directed in CSR) = %zu\n", G.n, G.col_idx.size());
    fprintf(stderr,"Using %d Pthreads\n", num_threads);

    auto t3 = clk::now();
    size_t iters = 0;
    vector<u32> labels = compute_connected_components(G, iters);
    auto t4 = clk::now();
    double cc_time = chrono::duration<double>(t4-t3).count();

    unordered_set<u32> comps(labels.begin(), labels.end());
    printf("Number of connected components: %zu\n", comps.size());
    printf("Number of iterations: %zu\n", iters);

    printf("Membership vector (label[v] for each vertex):\n");
    size_t limit = min((size_t)1000, labels.size());
    for(size_t i=0;i<limit;i++) printf("v%zu -> %u\n", i, labels[i]);

    printf("Timing info (seconds): I/O=%.3f, CSR build=%.3f, CC compute=%.3f\n",
           read_time, build_time, cc_time);

    return 0;
}
