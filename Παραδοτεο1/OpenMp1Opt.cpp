#include <bits/stdc++.h>
#include <chrono>
#include <omp.h>
using namespace std;
using u32 = uint32_t;
using u64 = uint64_t;
using clk = chrono::steady_clock;
//for friendster n mawi
// Edge in COO format (small temporary chunk)
struct MMEntry { u32 r, c; };

// CSR graph
struct CSR {
    u32 n = 0;
    vector<u64> row_ptr; // size n+1
    vector<u32> col_idx;
};

// Read header and get number of vertices, edges, and symmetric flag,same as the sequential part
bool read_matrix_market_header(const string &path, bool &symmetric, u32 &max_index) {
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
    max_index = max((u32)rows, (u32)cols) - 1;

    fclose(f);
    return true;
}

// Build CSR from MM file in chunks
CSR build_csr_chunked(const string &path, u32 n, bool symmetric) {
    CSR G; G.n = n;
    vector<u64> deg(n,0);

    //pass 1: Count degrees in chunks!
    const size_t CHUNK = 1000000;
    FILE *f = fopen(path.c_str(),"r");
    char buf[4096];
    // Skip header
    do { if (!fgets(buf, sizeof(buf), f)) { fclose(f); return G; } } while (buf[0]=='%');
    while (!feof(f)) {
        size_t count=0;
        vector<MMEntry> chunk;
        while(count<CHUNK && fgets(buf,sizeof(buf),f)){
            if(buf[0]=='%') continue;
            unsigned long long r=0,c=0;
            if(sscanf(buf,"%llu %llu",&r,&c)>=2){
                chunk.push_back({(u32)(r-1),(u32)(c-1)});
                count++;
            }
        }
        // Degree counting
        #pragma omp parallel for //each thread runs this code and all threads will enter the loop only for their assigned chunk indices
        for(size_t i=0;i<chunk.size();i++){ //also each thread has its own i
            u32 u=chunk[i].r, v=chunk[i].c;
            if(u!=v){  //for symmetric ones
                #pragma omp atomic  //atomic incremation  and write back to prevent race conditions between threads
                deg[u]++;
                #pragma omp atomic
                deg[v]++;
            } else if(!symmetric){  //If u==v (self-loop) and the graph is not symmetric,you only increment the degree once.
                #pragma omp atomic
                deg[u]++;
            }
        }
    }
    fclose(f);

    //Build row_ptr
    G.row_ptr.assign(n+1,0);
    for(u32 i=0;i<n;i++) G.row_ptr[i+1]=G.row_ptr[i]+deg[i];
    u64 m = G.row_ptr[n];
    G.col_idx.assign(m,0);

    //Current insertion index per row
    vector<u64> cur(n,0);
    for(u32 i=0;i<n;i++) cur[i]=G.row_ptr[i];

    //pass 2: Fill col_idx in chunks
    f = fopen(path.c_str(),"r");
    do { if (!fgets(buf, sizeof(buf), f)) { fclose(f); return G; } } while (buf[0]=='%');

    while (!feof(f)) {
        size_t count=0;
        vector<MMEntry> chunk;
        while(count<CHUNK && fgets(buf,sizeof(buf),f)){
            if(buf[0]=='%') continue;
            unsigned long long r=0,c=0;
            if(sscanf(buf,"%llu %llu",&r,&c)>=2){
                chunk.push_back({(u32)(r-1),(u32)(c-1)});
                count++;
            }
        }

        //Insert edges
        #pragma omp parallel for //creates team of threads 
        for(size_t i=0;i<chunk.size();i++){
            u32 u=chunk[i].r, v=chunk[i].c;
            if(u==v && !symmetric){
                u64 pos; 
                #pragma omp atomic capture //atomic actions to prevent race conditions between thread,this is a single atomic read-modify-write.
                { pos=cur[u]; cur[u]++; } //pos=cur[u] , cur[u]++; increment
                G.col_idx[pos]=v;
            } else if(u!=v){
                u64 pos_u,pos_v;
                #pragma omp atomic capture //same idea but here,we insert the undirected edge twice:
                { pos_u=cur[u]; cur[u]++; }
                #pragma omp atomic capture
                { pos_v=cur[v]; cur[v]++; }
                G.col_idx[pos_u]=v;
                G.col_idx[pos_v]=u;
            }
        }
    }
    fclose(f);

    //sort & deduplicate each row 
    #pragma omp parallel for schedule(dynamic,100) //each thread gets chunks of 100
    for(u32 i=0;i<n;i++){
        u64 s=G.row_ptr[i], t=G.row_ptr[i+1];
        if(t<=s+1) continue;
        sort(G.col_idx.begin()+s,G.col_idx.begin()+t);
        u64 write=s;
        for(u64 p=s;p<t;p++) if(p==s || G.col_idx[p]!=G.col_idx[p-1]) G.col_idx[write++]=G.col_idx[p];
        G.row_ptr[i+1]=write;
    }

    // Compaction
    vector<u32> new_cols; new_cols.reserve(G.col_idx.size());
    vector<u64> new_row_ptr(n+1,0);
    new_row_ptr[0]=0;
    for(u32 i=0;i<n;i++){
        u64 s=G.row_ptr[i], t=G.row_ptr[i+1];
        for(u64 p=s;p<t;p++) new_cols.push_back(G.col_idx[p]);
        new_row_ptr[i+1]=new_cols.size();
    }
    G.col_idx.swap(new_cols);
    G.row_ptr.swap(new_row_ptr);

    return G;
}

// Parallel label propagation for connected components
vector<u32> compute_connected_components(const CSR &G, size_t &num_iters){
    u32 n = G.n;
    vector<u32> labels(n);

    #pragma omp parallel for  //each thread writes to a different labels[i],so no races occur.
    for(u32 i=0;i<n;i++) labels[i]=i;

    bool changed=true;
    num_iters=0;

    while(changed){
        changed=false;
        num_iters++;

        #pragma omp parallel for schedule(dynamic,100) //each thread processes a subset of nodes
        for(u32 u=0;u<n;u++){ //each thread reads all neighbor labels of vertex u,computes the minimum label among them and if the minimum label is smaller than the current label, it updates labels[u]
            u32 min_label=labels[u];
            for(u64 p=G.row_ptr[u];p<G.row_ptr[u+1];p++){
                u32 v=G.col_idx[p];
                min_label=min(min_label, labels[v]);
            }
            if(min_label<labels[u]){
                labels[u]=min_label;
                #pragma omp atomic write //atomic because many threads may discover a change simultaneously,so we are safe from race conditions
                changed=true;
            }
        }
    }

    // Relabel to consecutive IDs
    vector<u32> unique_labels(labels.begin(),labels.end());
    sort(unique_labels.begin(),unique_labels.end());
    unique_labels.erase(unique(unique_labels.begin(),unique_labels.end()),unique_labels.end());
    unordered_map<u32,u32> label_map;
    for(u32 i=0;i<unique_labels.size();i++) label_map[unique_labels[i]]=i;
    #pragma omp parallel for //the same,team of threads for the mapping of our labels
    for(u32 i=0;i<n;i++) labels[i]=label_map[labels[i]];

    return labels;
}

int main(int argc, char **argv){
    if(argc<2){ fprintf(stderr,"Usage: %s graph.mtx [num_threads]\n",argv[0]); return 1; }
    string path=argv[1];
    if(argc>2) omp_set_num_threads(atoi(argv[2]));

    bool symmetric=false;
    u32 max_index=0;
    if(!read_matrix_market_header(path,symmetric,max_index)){ fprintf(stderr,"Failed to read header\n"); return 1; }
    u32 n = max_index+1;

    auto t0=clk::now();
    CSR G = build_csr_chunked(path,n,symmetric);
    auto t1=clk::now();
    double build_time=chrono::duration<double>(t1-t0).count();

    fprintf(stderr,"Graph: n=%u, edges=%zu\n",G.n,G.col_idx.size());
    fprintf(stderr,"Using %d OpenMP threads\n",omp_get_max_threads());
    fprintf(stderr,"CSR build time: %.3f s\n",build_time);

    auto t2=clk::now();
    size_t iters=0;
    vector<u32> labels = compute_connected_components(G,iters);
    auto t3=clk::now();
    double cc_time=chrono::duration<double>(t3-t2).count();

    unordered_set<u32> comps(labels.begin(),labels.end());
    printf("Number of connected components: %zu\n",comps.size());
    printf("Number of iterations: %zu\n",iters);

    // Print first 1000 nodes
    size_t limit=min((size_t)1000,labels.size());
    for(size_t i=0;i<limit;i++) printf("v%zu -> %u\n",i,labels[i]);

    printf("Timing (s): CSR build=%.3f, CC compute=%.3f\n",build_time,cc_time);

    return 0;
}
