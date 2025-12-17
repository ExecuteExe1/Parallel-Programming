#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <chrono>

using namespace std;
//using this to count time
using clk = chrono::steady_clock;

//constructor of CSRGraph
struct CSRGraph {
    int n; // vertices
    vector<uint64_t> rowptr; //row pointers: rowptr[i] = start of neighbors of vertex i in colind
    vector<int> colind;  //column indices: neighbors of all vertices stored consecutively
};

// Load Matrix Market graph in parallel 
CSRGraph load_mtx_as_csr(const string &filename) {
    ifstream fin(filename);//open MTX file
    if (!fin) { cerr << "Cannot open " << filename << endl; exit(1); }
    
    //Skip comments in Matrix Market file
    string line;
    while (getline(fin, line)) if (!line.empty() && line[0] != '%') break;

      //Read matrix dimensions: nrows = number of vertices
    int nrows, ncols, nnz;
    stringstream ss(line);
    ss >> nrows >> ncols >> nnz;
    
     //Build adjacency list
    vector<vector<int>> adj(nrows);

    int u, v;
    while (fin >> u >> v) {
        u--; v--; // 0-based
        adj[u].push_back(v);
        if (u != v) adj[v].push_back(u); //undirected graph
    }
    
  //Convert adjacency list to CSR format
    CSRGraph g;
    g.n = nrows;
    g.rowptr.resize(g.n + 1);
    g.rowptr[0] = 0;
    for (int i = 0; i < g.n; i++) g.rowptr[i+1] = g.rowptr[i] + adj[i].size();

    g.colind.resize(g.rowptr[g.n]);
    #pragma omp parallel for
    for (int i = 0; i < g.n; i++) {
        int idx = g.rowptr[i];
        for (int x : adj[i])
            g.colind[idx++] = x;
    }
    return g;
}

//Connected Components with hybrid MPI + OpenMP 
void cc_mpi_omp(const CSRGraph &G, vector<int> &label, int rank, int size) {
    int n = G.n;  //number of vertices
    label.resize(n);    //output labels (component ID) for each vertex
    vector<int> new_label(n); //temporary array to store updated labels

    //Initialize each vertex with its own label (label = vertex index)
    #pragma omp parallel for
    for (int v = 0; v < n; v++) label[v] = v;

    bool global_changed = true; //flag to check if any label changed globally

    while (global_changed) {
        bool local_changed = false; //flag to check if any label changed on this rank
        
         //Update labels for vertices assigned to this MPI rank
         //Each rank processes vertices: v = rank, rank + size, rank + 2*size, ...
         //OpenMP parallelizes this loop across threads

        #pragma omp parallel for reduction(||:local_changed)
        for (int v = rank; v < n; v += size) { // simple static partition
            int best = label[v]; //current label of vertex v
            
             //Iterate over neighbors of vertex v in CSR format
            for (uint64_t e = G.rowptr[v]; e < G.rowptr[v+1]; e++) {
                int u = G.colind[e];  //neighbor vertex
                if (label[u] < best) best = label[u]; //adopt minimum label among neighbors
            }
            new_label[v] = best;  //store updated label
            if (best != label[v]) local_changed = true; //mark change if label updated
        }

         //Copy labels for vertices not assigned to this rank
         //Ensures that unprocessed vertices retain their previous label
        #pragma omp parallel for
        for (int v = 0; v < n; v++)
            if (v % size != rank) new_label[v] = label[v];

        label.swap(new_label);  //efficiently update labels for next iteration

        // Global check across MPI ranks  
        MPI_Allreduce(&local_changed, &global_changed, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);
    }
}

//Main 
int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);   //Initialize the MPI environment

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //Get current process rank
    MPI_Comm_size(MPI_COMM_WORLD, &size); //Get total number of MPI processes

    if (argc < 2) {
        if (rank == 0) cerr << "Usage: mpirun -np N ./program graph.mtx" << endl;
        MPI_Finalize();   //Clean up MPI
        return 1;
    }

    string filename = argv[1];  //Graph file in Matrix Market format


    CSRGraph G;
    double t0 = MPI_Wtime();  //Start timing
    if (rank == 0) {         //Only rank 0 loads the graph from file
        cout << "Loading graph..." << endl; 
        G = load_mtx_as_csr(filename);   //Load MTX file and convert to CSR
        cout << "Graph: n=" << G.n << ", edges=" << G.colind.size() << endl;
    }

    // Broadcast graph size
    MPI_Bcast(&G.n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Resize on other ranks
    if (rank != 0) G.rowptr.resize(G.n + 1);
    MPI_Bcast(G.rowptr.data(), G.n + 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);

    uint64_t m = G.rowptr[G.n];
    if (rank != 0) G.colind.resize(m); 
    MPI_Bcast(G.colind.data(), m, MPI_INT, 0, MPI_COMM_WORLD);

    double t1 = MPI_Wtime();
    if (rank == 0) cout << "CSR loaded in " << t1 - t0 << " s" << endl;

    //Compute CC
    vector<int> label;   //Output labels for each vertex
    double t2 = MPI_Wtime(); //Start timing CC computation
    cc_mpi_omp(G, label, rank, size); //Hybrid MPI + OpenMP CC computation
    double t3 = MPI_Wtime();   //End timing

    if (rank == 0) {   
        unordered_set<int> comps(label.begin(), label.end());   //Count unique labels
        cout << "Connected components: " << comps.size() << endl;
        cout << "CC computation time: " << t3 - t2 << " s" << endl;
        cout << "Total execution time: " << t3 - t0 << " s" << endl;
    }

    // Optional: debug first 10 nodes
    for (int i = 0; i < min(10, (int)label.size()); i++)
        if (rank == 0) cout << "v" << i << " -> " << label[i] << endl;

    MPI_Finalize(); //Finalize the MPI environment
    return 0;
}
