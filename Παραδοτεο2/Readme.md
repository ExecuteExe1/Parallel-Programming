# GPU-Accelerated Graph Processing – Connected Components Detection

## Project Overview

This project was developed as part of the **Parallel Programming** course at the School of Electrical and Computer Engineering, Aristotle University of Thessaloniki (AUTH).

The objective is to accelerate **Connected Components (CC) detection in large-scale graphs using GPU parallel computing**, implemented in Julia with CUDA and KernelAbstractions.jl.

The project focuses on optimizing graph preprocessing, memory efficiency, GPU execution, and data transfer overhead, comparing performance against previous CPU-based implementations.

## Implementation

The application processes large graphs stored in **Matrix Market (.mtx)** format and converts them into a GPU-friendly **Compressed Sparse Row (CSR)** representation.

The implementation consists of four main stages:

1. **Graph Loading:** Efficiently parses Matrix Market files using Int32 representations, processes graph edges, and handles undirected graph connectivity.
2. **CSR Construction:** Builds compact adjacency structures with contiguous memory layouts suitable for efficient GPU memory access.
3. **Binary CSR Storage:** Saves preprocessed graph structures in binary format, eliminating expensive parsing and reconstruction during subsequent executions.
4. **GPU Connected Components Detection:** Uses CUDA kernels to perform iterative label propagation, assigning each vertex the minimum component label among itself and its neighbors until convergence.

## GPU Parallelization

The Connected Components algorithm is implemented using **Julia, CUDA.jl, and KernelAbstractions.jl**.

Each GPU thread processes a graph vertex, traversing its neighbors through CSR arrays and updating component labels.

A convergence flag tracks whether labels change between iterations, allowing the algorithm to terminate once the component assignments stabilize.

Additional experiments investigate different GPU thread configurations and their impact on execution time, memory bandwidth utilization, and overall performance.

## Experimental Evaluation

Performance was evaluated on three large-scale graph datasets:

- **com-Friendster**
- **com-Orkut**
- **com-LiveJournal**

The experiments measure graph loading time, GPU computation time, and total execution time.

### Performance Highlights

With binary CSR files already generated, the report presents the following results:

| Dataset | Total Execution Time | Speedup vs. Assignment 1 |
|---|---:|---:|
| com-Friendster | 1.476 s | ~3.8× |
| com-Orkut | 0.836 s | ~7× |
| com-LiveJournal | 1.083 s | ~7.95× |

The experiments demonstrate significant improvements through GPU acceleration and optimized graph storage.

These reported speedups compare against the best CPU implementations from Assignment 1 and include the benefit of reusing preprocessed binary CSR data.

## Performance Analysis

The project also investigates the relationship between GPU parallelism and memory bandwidth.

Experimental results show that increasing the number of GPU threads does not necessarily reduce execution time, particularly for memory-bound graph algorithms.

The analysis highlights the importance of:

- Efficient memory access patterns.
- Compact graph representations.
- Minimizing repeated preprocessing.
- Choosing appropriate GPU execution configurations.
- Understanding memory-bandwidth limitations.

## Technologies

**Julia · CUDA.jl · KernelAbstractions.jl · GPU Computing · CSR · Graph Algorithms · Connected Components · Parallel Programming · Performance Optimization**

## Conclusion

This project demonstrates how GPU parallel computing and memory-efficient graph representations can significantly accelerate large-scale graph processing.

By combining custom CUDA kernels, optimized CSR construction, and binary graph caching, the implementation achieves substantial reductions in execution time compared with earlier CPU-based approaches.

The experimental evaluation also illustrates the importance of memory bandwidth, data layout, and GPU thread configuration when optimizing graph algorithms.
