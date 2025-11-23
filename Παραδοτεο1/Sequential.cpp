#include <bits/stdc++.h>
#include <chrono>

using namespace std;
using u32 = uint32_t;
using u64 = uint64_t;
using clk = chrono::steady_clock;

//Struct for the COO
struct MMEntry { u32 r, c; };

//Reading the Matrix Market File
bool read_matrix_market_coo(const string &path, vector<MMEntry> &coords,  //file path to load,vector to store the rows/colums{edge pairs},boolean flag for symmetric or not,output for the largest index
                            bool &symmetric, u32 &max_index) {
    FILE *f = fopen(path.c_str(), "r"); //opens the file
    if (!f) { fprintf(stderr, "Cannot open '%s'\n", path.c_str()); return false; } //debugging line if there is any issue opening the file

    char buf[4096];  //huge buffer to read our file's lines
    if (!fgets(buf, sizeof(buf), f)) { fclose(f); return false; }  //if we fail to read the first line,the file closes immediately
    string header(buf); //converts string from c--->c++
    if (header.find("MatrixMarket") == string::npos) {  //checks if the header says that it is a Matrix Market file
        fprintf(stderr, "Not a Matrix Market file\n"); fclose(f); return false;  //if false we do print error message and close!
    }

    symmetric = (header.find("symmetric") != string::npos);  //if symmetric is detected flag=true otherwise false

    // Skip comments
    do { if (!fgets(buf, sizeof(buf), f)) { fclose(f); return false; } }
    while (buf[0] == '%');

    unsigned long long rows=0, cols=0, nnz=0; //variables to save our matrix's dimensions
    if (sscanf(buf, "%llu %llu %llu", &rows, &cols, &nnz) < 2) {   //Parses the numbers from the line scanned just before.if Parsing fails=error
        fprintf(stderr, "Bad size line: %s\n", buf); fclose(f); return false;
    }

    coords.clear(); //clears coordinate vector in case it had previous data for some unknown reason
    max_index = 0; 

    while (fgets(buf, sizeof(buf), f)) {  //while to read every line till the end!
        if (buf[0] == '%') continue;  //skip lines that begin with % just in case we have comments
        unsigned long long r=0,c=0;   //variables to store parsed rows and colums
        if (sscanf(buf, "%llu %llu", &r, &c) >= 2) {   //reads  the 2 integers from the line aka the row+colum
            if(r<1 || c<1) continue;   //if it is using 1 base indexing....skip! plus ignore negative or 0 entries
            coords.push_back({(u32)(r-1),(u32)(c-1)});  //Converts from  1-based → 0-based indexing and stores the edge
            max_index = max(max_index, max((u32)(r-1),(u32)(c-1)));  //Keeps track of the largest index seen so far
        }
    }

    fclose(f); //close the file

    if(symmetric){  
        size_t original_size = coords.size();  //Store how many entries existed originally (so we don’t reprocess new ones).
        for(size_t i=0;i<original_size;i++){  //for every original entry:
            if(coords[i].r != coords[i].c)  //if it is not a diagonal element (u != v):
                coords.push_back({coords[i].c, coords[i].r}); //add reverse edge (v,u)
        }
    }

    fprintf(stderr, "Read %zu entries. symmetric=%d, max_index=%u\n", //# of coordinate entries,symmetric or not,largest node index
            coords.size(), (int)symmetric, max_index);
    return true;
}

// CSR structure n=#of vertices in the graph
struct CSR {
    u32 n = 0;
    vector<u64> row_ptr; //prefix-sum array, size.Aka row_ptr[i] is the starting index in col_idx for neighbors of vertex i and row_ptr[i+1]=degree of vertex i 
                         //f.e row_ptr= 0          3      6      9
                         //    col_idx=[2,7,9] [0 4 5]  [1 6]  
                         //            ^neighbors of 0
    vector<u32> col_idx; //adjacency lists of all vertices stored back-to-back in one array.
};

// Sequential CSR builder
CSR build_csr_undirected(u32 n, vector<MMEntry> &coords) {
    vector<u64> deg(n,0);    //creates a degree array and inti to 0 it also counts how many edges touch each vertex

    // Sequential degree counting
    for (size_t i = 0; i < coords.size(); i++) {  // for each r,c we do
        auto &e = coords[i];  
        if(e.r != e.c) {  //skip if r==c
            deg[e.r]++;
            deg[e.c]++; //increase degree of both endpoints because the graph is undirected and each undirected edge apears TWICE  in adj lists r-->c and c-->r
        }
    }

    CSR G; G.n = n;  //we construct the row_ptr,allocate a CSR structure, set node count
    G.row_ptr.assign(n+1,0); //create row_ptr with size n+1 and init all to zeros
    for(u32 i=0;i<n;i++) G.row_ptr[i+1] = G.row_ptr[i] + deg[i]; //Builds a running prefix sum

    u64 m = G.row_ptr[n];  //m= total number of adjacency entries (twice the number of undirected edges)
    G.col_idx.assign(m,0); //big enough to hold all neighbor entries

    // Current position for each row
    vector<u64> cur = G.row_ptr; 

    // Sequential edge insertion,basically we loop over all original COO edges
    for(size_t idx=0; idx<coords.size(); idx++){
        auto &e = coords[idx];
        u32 u=e.r, v=e.c;
        if(u==v) continue; //skip self-loops

        u64 pos_u = cur[u]++;  //each edge adds two adjacency entries
        u64 pos_v = cur[v]++;  // curl basically points where to write next

        G.col_idx[pos_u] = v;
        G.col_idx[pos_v] = u; //Store neighbors into adjacency list
    }  //this step builds an unsorted CSR

    // Sequential sort + dedup
    for(u32 i=0;i<n;i++){  
        u64 s=G.row_ptr[i], t=G.row_ptr[i+1];
        if(t<=s+1) continue;//skips entries of <1

        sort(G.col_idx.begin()+s, G.col_idx.begin()+t);

        u64 write = s;
        for(u64 p=s;p<t;p++)   //classic sorted-array deduplication.
            if(p==s || G.col_idx[p]!=G.col_idx[p-1])
                G.col_idx[write++] = G.col_idx[p]; //writes unique elements back starting at write
                                                //write ends at new end of row
        G.row_ptr[i+1] = write; //shrinks row to remove duplicates
    }

    // Compact (sequential)   //after deduplication, some rows are shorter,we must pack everything tightly again.
    vector<u32> new_cols;
    new_cols.reserve(G.col_idx.size());  //temporary vector for compact adjacency list.

    vector<u64> new_row_ptr(n+1,0);  //new row_ptr will reflect compacted layout
    for(u32 i=0;i<n;i++){  //copy all unique adjacency entries sequentially
        u64 s = G.row_ptr[i], t = G.row_ptr[i+1];
        for(u64 p=s;p<t;p++) new_cols.push_back(G.col_idx[p]);
        new_row_ptr[i+1] = new_cols.size();
    }

    G.col_idx.swap(new_cols);
    G.row_ptr.swap(new_row_ptr);  //efficiently replace arrays without copying

    return G;
}

// Sequential CC via label propagation
vector<u32> compute_connected_components(const CSR &G, size_t &num_iters){ 
    u32 n = G.n;   //store number of vertices for convenience
    vector<u32> labels(n); //allocate labels[i] for all nodes.
    for(u32 i=0; i<n; i++) labels[i] = i;  //initialize each node’s label to itself,this means each node starts as its own component.
    num_iters = 0;
    bool changed = true;  //flag to detect convergence aka  labels stop changing

    while(changed) { //Main iteration loop aka keep going until no labels change
        changed = false;
        num_iters++;

        for(u32 u=0; u<n; u++) { //loop over all vertices
            u64 s = G.row_ptr[u], t = G.row_ptr[u+1]; //range of neighbors for node u
            u32 min_label = labels[u];  //start with the node's current label

            for(u64 p=s; p<t; p++) //find the minimum label among neighbors
                min_label = min(min_label, labels[G.col_idx[p]]);

            if(min_label < labels[u]) { //if any neighbor has a smaller label,update this node’s label
                labels[u] = min_label;
                changed = true;
            }
        }

        for(u32 u=0; u<n; u++) {                //A second "smoothing" pass that pushes labels even faster,this helps convergence but is still sequential
            u64 s = G.row_ptr[u], t = G.row_ptr[u+1];
            for(u64 p=s; p<t; p++)
                labels[u] = min(labels[u], labels[G.col_idx[p]]);
        }
    }

    // Normalize labels
    vector<u32> uniq = labels; //copy all labels.
    sort(uniq.begin(), uniq.end()); //
    uniq.erase(unique(uniq.begin(), uniq.end()), uniq.end()); //remove duplicates--> now only unique component IDs

    unordered_map<u32,u32> mp;  //map old labels
    for(u32 i=0;i<uniq.size();i++) mp[uniq[i]] = i; //assign new IDs

    for(u32 i=0;i<n;i++) labels[i] = mp[labels[i]];//rewrite each vertex’s label using compressed IDs.

    return labels;
}

// Sequential BFS CC
vector<u32> compute_connected_components_bfs(const CSR &G, size_t &num_iters){  //second CC implementation using BFS from each unvisited node
    u32 n = G.n; //#of vertices
    vector<u32> labels(n, UINT32_MAX);  //initialize all labels to "unvisited"
    num_iters = 0;
    u32 comp = 0; //current component number

    for(u32 v=0; v<n; v++){
        if(labels[v] != UINT32_MAX) continue; //skip already-visited nodes

        queue<u32> q; //start BFS from node v, assign it the component label
        q.push(v);
        labels[v] = comp;

        while(!q.empty()){ //pop the next node
            u32 u = q.front(); q.pop(); 
            u64 s=G.row_ptr[u], t=G.row_ptr[u+1]; //range of its neighbors

            for(u64 p=s;p<t;p++){  //visit each neighbor and mark it
                u32 nbr = G.col_idx[p];
                if(labels[nbr] == UINT32_MAX){
                    labels[nbr] = comp;
                    q.push(nbr);
                }
            }
            num_iters++; //increase iteration count per BFS pop.
        }
        comp++; //finished exploring the component — move to the next one
    }
    return labels; //return BFS component labels
}

int main(int argc, char **argv){
    if(argc<2){
        fprintf(stderr,"Usage: %s graph.mtx\n", argv[0]); //require input file path
        return 1;
    }

    string path=argv[1]; //store path

    //prepare variables for reading MatrixMarket file
    vector<MMEntry> coords;
    bool symmetric=false;
    u32 max_index=0;
    //start timing I/O
    auto t0 = clk::now();
    //load file as COO coordinate list
    if(!read_matrix_market_coo(path, coords, symmetric, max_index)){
        fprintf(stderr,"Failed to read file\n");
        return 1;
    }

    u32 n = max_index + 1;
    //end I/o timing
    auto t1 = clk::now();
    //convert COO → CSR
    CSR G = build_csr_undirected(n, coords);
    //start the CSR timing
    auto t2 = clk::now();
   //compute timings
    double read_time = chrono::duration<double>(t1-t0).count();
    double build_time = chrono::duration<double>(t2-t1).count();
    
    fprintf(stderr,"I/O time: %.3f s, CSR build: %.3f s\n", read_time, build_time);
    fprintf(stderr,"Graph: n=%u, edges stored=%zu\n", G.n, G.col_idx.size());

     
    size_t iters=0;
    auto t3 = clk::now();
    vector<u32> labels = compute_connected_components(G, iters);
    auto t4 = clk::now();
    //compute CC time
    double cc_time = chrono::duration<double>(t4-t3).count();
  
    unordered_set<u32> comps(labels.begin(), labels.end());
    printf("Connected components: %zu\n", comps.size());
    printf("Iterations: %zu\n", iters);

    printf("Membership (first 1000 nodes):\n");
    for(size_t i=0;i<min<size_t>(1000, labels.size());i++)
        printf("v%zu -> %u\n", i, labels[i]);

    printf("Timing (seconds): I/O=%.3f, Build=%.3f, CC=%.3f\n",
           read_time, build_time, cc_time);

    return 0;
}
