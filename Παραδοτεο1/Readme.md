# Parallel Graph Processing – Connected Components Detection

## Project Overview

This project was developed as part of the **Parallel Programming** course at the School of Electrical and Computer Engineering, Aristotle University of Thessaloniki (AUTH).

The objective is to implement and optimize **Connected Components (CC) detection in large-scale graphs**, comparing sequential execution against three parallel programming frameworks:

- **OpenMP**
- **OpenCilk**
- **POSIX Threads (Pthreads)**

The project investigates parallelization strategies, memory management, workload distribution, and execution-time improvements across graph datasets with different sizes and connectivity characteristics.

## Implementation

The application reads graphs stored in **Matrix Market (.mtx)** format and constructs a memory-efficient **Compressed Sparse Row (CSR)** representation.

Connected components are identified using **Breadth-First Search (BFS)**, assigning a component identifier to vertices belonging to the same connected component.

The implementation consists of three main stages:

1. **Graph Loading:** Reading Matrix Market files, handling symmetric edges, and preparing the graph data.
2. **CSR Construction:** Building a compact adjacency representation, sorting neighbors, and removing duplicate edges.
3. **Connected Components Detection:** Traversing the graph using BFS to identify and label disconnected components.

## Parallel Implementations

### OpenMP

Parallelizes graph preprocessing, CSR construction, and selected connected-components operations.

Different implementations were developed to accommodate graph size and memory requirements, including chunk-based processing for memory-intensive datasets.

### OpenCilk

Uses task-based parallel programming to optimize graph processing.

Alternative strategies were explored, including parallel neighbor sorting and selectively retaining sequential connected-components traversal when parallelization introduced unnecessary overhead.

### POSIX Threads (Pthreads)

Implements explicit multithreading for CSR construction, including parallel vertex-degree counting and neighbor insertion using atomic operations.

The connected-components stage is kept sequential in the evaluated implementation, based on observations that parallel traversal did not consistently improve performance.

## Experimental Evaluation

Performance was evaluated using four graph datasets:

- com-Orkut
- com-LiveJournal
- com-Friendster
- MAWI

Each implementation was benchmarked against a sequential baseline, measuring execution time for graph loading, CSR construction, connected-components detection, and total processing.

### Performance Highlights

The experiments demonstrated workload-dependent improvements, with reported speedups including:

| Framework | Dataset | Speedup |
|---|---|---|
| OpenMP | com-Orkut | 1.43× |
| OpenMP | com-LiveJournal | 1.65× |
| OpenCilk | com-LiveJournal | 1.88× |
| OpenCilk | MAWI | 1.67× |
| Pthreads | com-LiveJournal | 1.446× |
| Pthreads | MAWI | ~1.8× |

Results also demonstrated that **parallelization does not always guarantee better performance**. Thread synchronization, memory overhead, graph structure, and workload characteristics can introduce additional costs, making hybrid sequential-parallel approaches preferable in certain cases.

## Technologies

**C++ · OpenMP · OpenCilk · POSIX Threads · Graph Algorithms · BFS · CSR · Parallel Computing · Performance Optimization**

## Conclusion

This project demonstrates the practical challenges of parallel graph processing and highlights the importance of selecting appropriate parallelization strategies based on workload characteristics.

The experimental comparison provides insights into thread-level parallelism, memory efficiency, scalability, and the trade-offs between parallel execution overhead and computational performance.
