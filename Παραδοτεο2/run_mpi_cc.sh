#!/bin/bash
#SBATCH --partition=rome
#SBATCH --time=00:10:00
#SBATCH --job-name=CILK_MPI_test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=60
#SBATCH --mem=64GB

module load gcc openmpi

# Compile
mpicxx -O3 -fopenmp cc_c_parallel.cpp -o cc_c_parallel

# Force MPI/TCP to avoid UCX timeouts
export UCX_TLS=tcp
export UCX_WARN_UNUSED_ENV_VARS=n
export OMPI_MCA_btl=tcp,self

# Set OpenMP threads
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Graph file
if [ $# -eq 0 ]; then
    echo "No graph file specified, defaulting to test.mtx"
    GRAPH_FILE="test.mtx"
else
    GRAPH_FILE="$1"
fi

if [ ! -f "$GRAPH_FILE" ]; then
    echo "Error: Graph file '$GRAPH_FILE' not found!"
    exit 1
fi

echo "Running with graph: $GRAPH_FILE"
echo "MPI tasks: $SLURM_NTASKS, OpenMP threads per task: $OMP_NUM_THREADS"

# Run
srun --ntasks=$SLURM_NTASKS ./cc_c_parallel "$GRAPH_FILE"
