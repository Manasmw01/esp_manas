// Copyright (c) 2011-2024 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0

#ifndef __DATATYPES__
#define __DATATYPES__

#include "ac_int.h"
#include "ac_fixed.h"
#include "kalman_filter_specs.hpp"

#define N const_mat_dim
#define FPDATA_WL DATA_WIDTH
// #define FPDATA_IL DATA_WIDTH/2
#define FPDATA_IL 6

typedef ac_int<DMA_WIDTH> DMA_WORD;
typedef ac_int<FPDATA_WL> FPDATA_WORD;
typedef ac_fixed<FPDATA_WL, FPDATA_IL> FPDATA;

inline void int2fx(const FPDATA_WORD& in, FPDATA& out)
{ out.set_slc(0,in.slc<FPDATA_WL>(0)); }

inline void fx2int(const FPDATA& in, FPDATA_WORD& out)
{ out.set_slc(0,in.slc<FPDATA_WL>(0)); }




#ifdef PRINT_STATEMENTS
inline void print_matrix(FPDATA matrixx[N][N], uint32_t kalman_mat_dim)
{
    for (int i = 0; i < kalman_mat_dim; i++)
        for (int j = 0; j < kalman_mat_dim; j++)
            std::cout << std::setw(20) << matrixx[i][j] << ((j == kalman_mat_dim - 1) ? "\n": "\t");
    std::cout << "\n";
}



inline void print_vector(FPDATA vec[N] , uint32_t kalman_mat_dim)
{
    for (int i = 0; i < kalman_mat_dim; i++)
            std::cout << std::setw(20) << vec[i] << "\t";
    std::cout << "\n";
}
// Functions to print matrices and vectors
inline void print_matrix_new(FPDATA* matrix, int rows, int cols) {
    // printf("Matrix (%d x %d):\n", rows, cols);
    for (int i = 0; i < rows; i++) {
        // printf("(Row %d)\t:", i);
        for (int j = 0; j < cols; j++) {
              cout << matrix[i * cols + j] << "\t"; 
        }
        cout << std::endl;
    }
}
#endif

inline void copymat(FPDATA A[N][N], FPDATA result[N][N], uint32_t kalman_mat_dim) {
    for (int i = 0; i < kalman_mat_dim; i++) {
        for (int j = 0; j < kalman_mat_dim; j++) {
            result[i][j] = A[i][j];
        }
    }
}

// inline void matrix_multiply(FPDATA A[MEAS_SIZE], FPDATA B[MEAS_SIZE], FPDATA C[MEAS_SIZE], uint32_t n, uint32_t m, uint32_t p) {
inline void matrix_multiply(FPDATA* A, FPDATA* B, FPDATA* C, uint32_t n, uint32_t m, uint32_t p) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < p; j++) {
            C[i * p + j] = 0;
            for (int k = 0; k < m; k++) {
                C[i * p + j] += A[i * m + k] * B[k * p + j];
            }
        }
    }
}



// Utility function implementations
inline void gauss_inverse(FPDATA* A, FPDATA* A_inv, int n) {
    // Augmenting the matrix A with identity matrix of same dimensions
    FPDATA augmented[n * 2 * n];

    // Create the augmented matrix
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            augmented[i * 2 * n + j] = A[i * n + j];  // A portion
            augmented[i * 2 * n + (j + n)] = (i == j) ? 1.0 : 0.0;  // Identity portion
        }
    }

    // Applying Gauss-Jordan elimination
    for (int i = 0; i < n; i++) {
        // Find the pivot row
        int pivot_row = i;
        for (int j = i + 1; j < n; j++) {
            if (augmented[j * 2 * n + i] > augmented[pivot_row * 2 * n + i]) {
                pivot_row = j;
            }
        }

        // Swap rows i and pivot_row
        if (pivot_row != i) {
            for (int k = 0; k < 2 * n; k++) {
                FPDATA temp = augmented[i * 2 * n + k];
                augmented[i * 2 * n + k] = augmented[pivot_row * 2 * n + k];
                augmented[pivot_row * 2 * n + k] = temp;
            }
        }

        // Make the diagonal elements 1
        FPDATA pivot = augmented[i * 2 * n + i];
        for (int k = 0; k < 2 * n; k++) {
            augmented[i * 2 * n + k] /= pivot;
        }

        // Make other elements in the column 0
        for (int j = 0; j < n; j++) {
            if (j != i) {
                FPDATA factor = augmented[j * 2 * n + i];
                for (int k = 0; k < 2 * n; k++) {
                    augmented[j * 2 * n + k] -= factor * augmented[i * 2 * n + k];
                }
            }
        }
    }

    // Copy the inverse from the augmented matrix
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A_inv[i * n + j] = augmented[i * 2 * n + (j + n)];
        }
    }
}


inline void matrix_add(FPDATA* A, FPDATA* B, FPDATA* C, int n, int m) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            C[i * m + j] = A[i * m + j] + B[i * m + j];
        }
    }
}

inline void matrix_transpose(FPDATA* A, FPDATA* AT, int n, int m) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            AT[j * n + i] = A[i * m + j];
        }
    }
}

inline void matrix_subtract(FPDATA* A, FPDATA* B, FPDATA* C, int n) {
    for (int i = 0; i < n; i++) {
        C[i] = A[i] - B[i];
    }
}

#endif
