#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>

#define min(a, b) ((a) < (b) ? (a) : (b))

// Quickselect-based function to find k smallest elements
void quickselect(double *distances, int *indices, int left, int right, int k);

// Function to calculate k nearest neighbors
void knnsearch(double *C, double *Q, int num_points, int dimensions, int k, int *indices, double *distances) {
    // Single-threaded calculation of distance matrix
    double *dist_matrix = (double *)malloc(num_points * num_points * sizeof(double));

    for (int j = 0; j < num_points; ++j) {
        for (int l = 0; l < num_points; ++l) {
            double dist = 0.0;
            // Compute squared distance between point j in Q and point l in C
            for (int d = 0; d < dimensions; ++d) {
                double diff = Q[j * dimensions + d] - C[l * dimensions + d];
                dist += diff * diff;
            }
            dist_matrix[j * num_points + l] = sqrt(dist); // Store the Euclidean distance
        }
    }

    // Find k smallest distances for each point in Q
    for (int j = 0; j < num_points; ++j) {
        double *dist_block = dist_matrix + j * num_points;
        int *idx_block = (int *)malloc(num_points * sizeof(int));
        for (int m = 0; m < num_points; ++m) idx_block[m] = m;

        quickselect(dist_block, idx_block, 0, num_points - 1, k);

        for (int m = 0; m < k; ++m) {
            indices[j * k + m] = idx_block[m];
            distances[j * k + m] = dist_block[idx_block[m]];
        }
        free(idx_block);
    }

    free(dist_matrix);
}

// Helper function for quickselect to partition array
int partition(double *distances, int *indices, int left, int right, int pivot_idx) {
    double pivot_val = distances[pivot_idx];
    int temp = indices[pivot_idx];
    indices[pivot_idx] = indices[right];
    indices[right] = temp;

    double temp_dist = distances[pivot_idx];
    distances[pivot_idx] = distances[right];
    distances[right] = temp_dist;

    int store_idx = left;
    for (int i = left; i < right; i++) {
        if (distances[i] < pivot_val) {
            temp = indices[i];
            indices[i] = indices[store_idx];
            indices[store_idx] = temp;

            temp_dist = distances[i];
            distances[i] = distances[store_idx];
            distances[store_idx] = temp_dist;

            store_idx++;
        }
    }
    temp = indices[store_idx];
    indices[store_idx] = indices[right];
    indices[right] = temp;

    temp_dist = distances[store_idx];
    distances[store_idx] = distances[right];
    distances[right] = temp_dist;

    return store_idx;
}

// Quickselect main function
void quickselect(double *distances, int *indices, int left, int right, int k) {
    if (left == right) return;

    int pivot_idx = left + rand() % (right - left + 1);
    pivot_idx = partition(distances, indices, left, right, pivot_idx);

    if (k == pivot_idx) return;
    else if (k < pivot_idx) quickselect(distances, indices, left, pivot_idx - 1, k);
    else quickselect(distances, indices, pivot_idx + 1, right, k);
}

int main() {
    // Example parameters
    int num_points = 10000;      // Number of points in C and Q
    int dimensions = 64;         // Dimensionality of each point
    int k = 10;                  // Number of nearest neighbors to find

    // Allocate memory for C and Q
    double *C = (double *)malloc(num_points * dimensions * sizeof(double));
    double *Q = (double *)malloc(num_points * dimensions * sizeof(double));

    // Example initialization (you may replace this with actual data)
    for (int i = 0; i < num_points * dimensions; ++i) {
        C[i] = rand() / (double)RAND_MAX;
        Q[i] = rand() / (double)RAND_MAX;
    }

    // Allocate memory for output arrays
    int *indices = (int *)malloc(num_points * k * sizeof(int)); // Output indices
    double *distances = (double *)malloc(num_points * k * sizeof(double)); // Output distances

    // Call the knnsearch function
    knnsearch(C, Q, num_points, dimensions, k, indices, distances);

    // Example output (print the first point's k-nearest neighbors)
    printf("Nearest neighbors for the first query point:\n");
    for (int i = 0; i < k; ++i) {
        printf("Index: %d, Distance: %f\n", indices[i], distances[i]);
    }

    // Free allocated memory
    free(C);
    free(Q);
    free(indices);
    free(distances);

    return 0;
}