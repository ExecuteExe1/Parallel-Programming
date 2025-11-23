#include <bits/stdc++.h>
#include <chrono>
#include <cilk/cilk.h>
#include <atomic>
using namespace std;
using u32 = uint32_t;
using u64 = uint64_t;
using clk = chrono::steady_clock;
//mawi and friendster

// CSR graph
struct CSR {
    u32 n;
    vector<u64> row_ptr;
    vector<u32> col_idx;
};

//Read Matrix Market header and get n, num_edges,sequential
bool read_matrix_market_detailed(const string &path, u32 &n, u64 &num_edges) {
    FILE *f = fopen(path.c_str(), "r");
    if (!f) { fprintf(stderr, "Cannot open '%s'\n", path.c_str()); return false; }

    char buf[4096];
    if (!fgets(buf, sizeof(buf), f)) { fclose(f); return false; }

    // Skip comments
    while (fgets(buf, sizeof(buf), f) && buf[0] == '%');

    unsigned long long rows=0, cols=0, nnz=0;
    if (sscanf(buf, "%llu %llu %llu", &rows, &cols, &nnz) < 2) {
        fprintf(stderr, "Bad size line: %s\n", buf); fclose(f); return false;
    }

    n = max(rows, cols);
    num_edges = nnz;
    fclose(f);
    fprintf(stderr, "Graph: n=%u, edges=%llu\n", n, nnz);
    return true;
}

//Build CSR with separate I/O and processing timing
CSR build_csr_with_timing(const string &path, u32 n, double &io_time, double &processing_time) {
    CSR G; G.n = n;
    vector<u32> out_deg(n, 0);

    char buf[4096];
    FILE *f = fopen(path.c_str(), "r");
    if (!f) { fprintf(stderr,"Cannot open file\n"); exit(1); }

    // Skip header and comments
    do { if (!fgets(buf, sizeof(buf), f)) break; } while (buf[0] == '%');

   
    //pass 1:count out-degrees (I/O)
    io_time = 0.0;
    auto io_start = clk::now();
    u64 line_count = 0, edges_found = 0;
    while (fgets(buf, sizeof(buf), f)) {
        if (buf[0]=='%') continue;
        unsigned long long r,c;
        if (sscanf(buf,"%llu %llu",&r,&c)>=2 && r>=1 && c>=1) {
            u32 u = r-1;
            if(u<n) { out_deg[u]++; edges_found++; }
        }
        line_count++;
        if(line_count%10000000==0) fprintf(stderr,"Pass1: processed %lu lines, %lu edges\n",line_count,edges_found);
    }
    auto io_end = clk::now();
    io_time += chrono::duration<double>(io_end-io_start).count();
    fprintf(stderr,"Pass1 complete: %lu edges counted\n", edges_found);

    
    //PROCESS: BUILD ROW_PTR
    auto proc_start = clk::now();
    G.row_ptr.resize(n+1);
    G.row_ptr[0]=0;
    for(u32 i=0;i<n;i++) G.row_ptr[i+1]=G.row_ptr[i]+out_deg[i];
    u64 total_edges = G.row_ptr[n];
    fprintf(stderr,"CSR row_ptr built: %lu edges\n", total_edges);
    processing_time = chrono::duration<double>(clk::now()-proc_start).count();

    //Allocate col_idx
    try { G.col_idx.resize(total_edges); } 
    catch(const std::bad_alloc &e) { 
        fprintf(stderr,"Memory allocation failed: %s\n", e.what()); exit(1); 
    }

    
    //pass 2: FILL COL_IDX (I/O)
    vector<u32> pos(n);
    for(u32 i=0;i<n;i++) pos[i]=G.row_ptr[i];

    fseek(f,0,SEEK_SET);
    do { if(!fgets(buf,sizeof(buf),f)) break; } while(buf[0]=='%');

    io_start = clk::now();
    line_count = 0; edges_found = 0;
    while (fgets(buf,sizeof(buf),f)) {
        if(buf[0]=='%') continue;
        unsigned long long r,c;
        if(sscanf(buf,"%llu %llu",&r,&c)>=2 && r>=1 && c>=1) {
            u32 u = r-1, v = c-1;
            if(u<n && v<n) { G.col_idx[pos[u]++]=v; edges_found++; }
        }
        line_count++;
        if(line_count%10000000==0) fprintf(stderr,"Pass2: processed %lu lines, %lu edges\n",line_count,edges_found);
    }
    io_end = clk::now();
    io_time += chrono::duration<double>(io_end-io_start).count();
    fclose(f);
    fprintf(stderr,"Pass2 complete: %lu edges inserted\n", edges_found);

    //PROCESS: SORT NEIGHBORS (processing)
    proc_start = clk::now();
    cilk_for(u32 i=0;i<n;i++) {
        u64 s=G.row_ptr[i], t=G.row_ptr[i+1];
        if(t>s) sort(G.col_idx.begin()+s,G.col_idx.begin()+t);
    }
    processing_time += chrono::duration<double>(clk::now()-proc_start).count();
    fprintf(stderr,"CSR neighbor sorting complete\n");

    return G;
}

//sequential BFS for weakly connected components
vector<u32> compute_cc_bfs(const CSR &G, size_t &num_iters) {
    u32 n = G.n;
    vector<u32> labels(n, UINT32_MAX);
    num_iters=0;
    for(u32 v=0;v<n;v++) {
        if(labels[v]!=UINT32_MAX) continue;
        queue<u32> q;
        q.push(v);
        labels[v]=v;
        while(!q.empty()) {
            u32 u=q.front(); q.pop();
            for(u64 p=G.row_ptr[u];p<G.row_ptr[u+1];p++){
                u32 nbr = G.col_idx[p];
                if(labels[nbr]==UINT32_MAX) { labels[nbr]=v; q.push(nbr); }
            }
            num_iters++;
        }
        if(v%1000000==0 && v>0) fprintf(stderr,"CC BFS: processed %u vertices\n",v);
    }
    return labels;
}

int main(int argc,char **argv){
    if(argc<2){ fprintf(stderr,"Usage: %s graph.mtx\n",argv[0]); return 1; }

    string path = argv[1];
    fprintf(stderr,"=== STARTING GRAPH PROCESSING ===\n");

    u32 n; u64 num_edges;
    auto t0 = clk::now();
    if(!read_matrix_market_detailed(path,n,num_edges)) return 1;
    auto t1 = clk::now();
    double header_time = chrono::duration<double>(t1-t0).count();

    double io_time=0, proc_time=0;
    auto t2 = clk::now();
    CSR G = build_csr_with_timing(path,n,io_time,proc_time);
    auto t3 = clk::now();
    double build_total = chrono::duration<double>(t3-t2).count();

    size_t iters=0;
    auto t4 = clk::now();
    vector<u32> labels = compute_cc_bfs(G,iters);
    auto t5 = clk::now();
    double cc_time = chrono::duration<double>(t5-t4).count();

    double total_time = chrono::duration<double>(t5-t0).count();

    // Count sampled components
    unordered_set<u32> comps;
    size_t sample = min((size_t)1000000, labels.size());
    for(size_t i=0;i<sample;i++) comps.insert(labels[i]);

    printf("\n=== FINAL RESULTS ===\n");
    printf("Graph: %u nodes, %zu edges\n", n, G.col_idx.size());
    printf("Connected components (sample first %zu vertices): %zu\n", sample, comps.size());
    printf("Iterations: %zu\n", iters);
    printf("\n=== TIMING ===\n");
    printf("Header reading:        %.3f s\n", header_time);
    printf("File I/O (both passes):%.3f s\n", io_time);
    printf("CSR processing:        %.3f s\n", proc_time);
    printf("CSR build total:       %.3f s\n", build_total);
    printf("CC computation:        %.3f s\n", cc_time);
    printf("TOTAL TIME:            %.3f s\n", total_time);

    return 0;
}
