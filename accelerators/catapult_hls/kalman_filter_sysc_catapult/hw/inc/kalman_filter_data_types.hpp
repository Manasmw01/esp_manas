// Copyright (c) 2011-2024 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0

#ifndef __DATATYPES__
#define __DATATYPES__

#include "ac_int.h"
#include "ac_fixed.h"
#include "kalman_filter_specs.hpp"
#include "ac_float.h"

#define N const_mat_dim
#define FPDATA_WL DATA_WIDTH

#define FPDATA_WL DATA_WIDTH
#define FPDATA_IL DATA_WIDTH/2


typedef ac_int<DMA_WIDTH> DMA_WORD;
typedef ac_int<FPDATA_WL> FPDATA_WORD;
typedef ac_fixed<FPDATA_WL, FPDATA_IL> FPDATA;


typedef ac_float<23, 0, 8> FLOAT_TYPE;

typedef FLOAT_TYPE FN_DATATYPE;

inline void int2fp(const FPDATA_WORD& in, FLOAT_TYPE& out) {
    float temp = *reinterpret_cast<const float*>(&in);  // reinterpret cast
    out = FLOAT_TYPE(temp);  // Assign float to FLOAT_TYPE
}

inline void fp2int(const FLOAT_TYPE& in, FPDATA_WORD& out) {
    float temp = in.to_float();  // Convert FLOAT_TYPE to float
    out = *reinterpret_cast<FPDATA_WORD*>(&temp);  // reinterpret cast to FPDATA_WORD
}


inline void int2fx(const FPDATA_WORD& in, FPDATA& out)
{ out.set_slc(0,in.slc<FPDATA_WL>(0)); }

inline void fx2int(const FPDATA& in, FPDATA_WORD& out)
{ out.set_slc(0,in.slc<FPDATA_WL>(0)); }




#ifdef PRINT_STATEMENTS
inline void print_matrix(FN_DATATYPE matrixx[N][N], uint32_t kalman_mat_dim)
{
    // for (int i = 0; i < kalman_mat_dim; i++)
    //     for (int j = 0; j < kalman_mat_dim; j++)
    //         std::cout << std::setw(20) << matrixx[i][j] << ((j == kalman_mat_dim - 1) ? "\n": "\t");
    // std::cout << "\n";
}



inline void print_vector(FN_DATATYPE vec[N] , uint32_t kalman_mat_dim)
{
    // for (int i = 0; i < kalman_mat_dim; i++)
    //         std::cout << std::setw(20) << vec[i] << "\t";
    // std::cout << "\n";
}
// Functions to print matrices and vectors
inline void print_matrix_new(FN_DATATYPE* matrix, int rows, int cols) {
    // printf("Matrix (%d x %d):\n", rows, cols);
    for (int i = 0; i < rows; i++) {
        // printf("(Row %d)\t:", i);
        for (int j = 0; j < cols; j++) {
            // cout << std::setprecision(30) << matrix[i * cols + j].to_float() << " "; 
            printf("%.30f ", matrix[i * cols + j].to_float());
        }
        cout << std::endl;
    }
}
#endif

inline void copymat(FN_DATATYPE A[N][N], FN_DATATYPE result[N][N], uint32_t kalman_mat_dim) {
    for (int i = 0; i < kalman_mat_dim; i++) {
        for (int j = 0; j < kalman_mat_dim; j++) {
            result[i][j] = A[i][j];
        }
    }
}

// inline void matrix_multiply(FPDATA A[MEAS_SIZE], FPDATA B[MEAS_SIZE], FPDATA C[MEAS_SIZE], uint32_t n, uint32_t m, uint32_t p) {
inline void matrix_multiply(FN_DATATYPE* A, FN_DATATYPE* B, FN_DATATYPE* C, uint32_t n, uint32_t m, uint32_t p) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < p; j++) {
            C[i * p + j] = 0;
            for (int k = 0; k < m; k++) {
                C[i * p + j] += A[i * m + k] * B[k * p + j];
            }
        }
    }
}


inline void inverse_clean(FN_DATATYPE new_mat[MEAS_SIZE][MEAS_SIZE], FN_DATATYPE out[MEAS_SIZE][MEAS_SIZE])
{

		 FN_DATATYPE ratio;
		 int i,j,k;

		 if(MEAS_SIZE == 2){
			 FN_DATATYPE a = new_mat[0][0];
			 FN_DATATYPE b = new_mat[0][1];
			 FN_DATATYPE c = new_mat[1][0];
			 FN_DATATYPE d = new_mat[1][1];

			 FN_DATATYPE det = FN_DATATYPE((a.to_float() * d.to_float()) - (b.to_float() * c.to_float()));

			 if (det == 0) {
			    //  printf("The matrix is not invertible.\n");
			     return;
			 }

			 out[0][0] = d / det;
			 out[0][1] = (-1) * b / det;
			 out[1][0] = (-1) * c / det;
			 out[1][1] = a / det;

			 return;
		 }

		 /* Applying Gauss Jordan Elimination */
		 for(i = 0; i < MEAS_SIZE; i++)
		 {
			  for(j = 0; j < MEAS_SIZE; j++)
			  {
				   if(i != j)
				   {
					    ratio = new_mat[j][i]/new_mat[i][i];
					    for(k = 0; k < MEAS_SIZE; k++)
					    {

					    	if(i == MEAS_SIZE-1){
					    		if(k == 0){//Calc the diagonal element first
					    			new_mat[j][j] = FN_DATATYPE(new_mat[j][j].to_float() - ratio.to_float()*new_mat[i][j].to_float());
					    			//out[j][j] = (out[j][j] - ratio*out[i][j]) / new_mat[j][j];
					    		}
					    		else if(k == j){
					    			new_mat[j][0] = FN_DATATYPE(new_mat[j][0].to_float() - ratio.to_float()*new_mat[i][0].to_float());
					    			//out[j][0] = (out[j][0] - ratio*out[i][0]) / new_mat[j][j];
					    		}
					    		else{
					    			new_mat[j][k] = FN_DATATYPE(new_mat[j][k].to_float() - ratio.to_float()*new_mat[i][k].to_float());
					    			//out[j][k] = (out[j][k] - ratio*out[i][k]) / new_mat[j][j];
					    		}

					    		out[j][k] = FN_DATATYPE((out[j][k].to_float() - ratio.to_float()*out[i][k].to_float()) / new_mat[j][j].to_float());
					    		//out[j][j] = out[j][j] / new_mat[j][j];
					    	}
					    	else{

								new_mat[j][k] = FN_DATATYPE(new_mat[j][k].to_float() - ratio.to_float()*new_mat[i][k].to_float());

								if(i > 0)
									out[j][k] = FN_DATATYPE(out[j][k].to_float() - ratio.to_float()*out[i][k].to_float());
								else{ //(i == 0)
									if(i == k)
										out[i][k] = 1;
									else
										out[i][k] = 0;
									if(j == k){
										if(i == k)
											out[j][k] = FN_DATATYPE(1 - ratio.to_float());
										else
											out[j][k] = 1;
										//out[j][k] = 1 - ratio*out[i][k];
									}
									else{
										if(i == k)
											out[j][k] = -ratio;
										else
											out[j][k] = 0;
										//out[j][k] = 0 - ratio*out[i][k];
									}
								}
					    	}

					    }

				   }

			  }
		 }
		 for(i = 0; i < MEAS_SIZE; i++)
		 {
			 out[MEAS_SIZE-1][i] = out[MEAS_SIZE-1][i] / new_mat[MEAS_SIZE-1][MEAS_SIZE-1];
		 }
		 return;
}

// Utility function implementations

inline void gauss_inverse(FN_DATATYPE* A, FN_DATATYPE* A_inv, int n) {
    // Augmenting the matrix A with identity matrix of same dimensions
    std::cout << "In Gauss Inverse\n";
    FN_DATATYPE augmented[N * 2 * N];

    // Create the augmented matrix
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            augmented[i * 2 * n + j] = A[i * n + j];  // A portion
            augmented[i * 2 * n + (j + n)] = (i == j) ? 1 : 0;  // Identity portion
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
                FN_DATATYPE temp = augmented[i * 2 * n + k];
                augmented[i * 2 * n + k] = augmented[pivot_row * 2 * n + k];
                augmented[pivot_row * 2 * n + k] = temp;
            }
        }

        // Make the diagonal elements 1
        FN_DATATYPE pivot = augmented[i * 2 * n + i];
        for (int k = 0; k < 2 * n; k++) {
            augmented[i * 2 * n + k] /= pivot;
        }

        // Make other elements in the column 0
        for (int j = 0; j < n; j++) {
            if (j != i) {
                FN_DATATYPE factor = augmented[j * 2 * n + i];
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


inline void matrix_add(FN_DATATYPE* A, FN_DATATYPE* B, FN_DATATYPE* C, int n, int m) {
// inline void matrix_add(FN_DATATYPE A[MEAS_SIZE][MEAS_SIZE], FN_DATATYPE B[MEAS_SIZE][MEAS_SIZE], FN_DATATYPE C[MEAS_SIZE][MEAS_SIZE], int n, int m) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            // C[i * m + j] = A[i * m + j] + B[i * m + j];
            C[i * m + j] = FLOAT_TYPE(A[i * m + j].to_float() + B[i * m + j].to_float());
        }
    }
}

inline void matrix_transpose(FN_DATATYPE* A, FN_DATATYPE* AT, int n, int m) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            AT[j * n + i] = A[i * m + j];
        }
    }
}

inline void matrix_subtract(FN_DATATYPE* A, FN_DATATYPE* B, FN_DATATYPE* C, int n) {
    for (int i = 0; i < n; i++) {
        C[i] = FLOAT_TYPE(A[i].to_float() - B[i].to_float());
    }
}

#endif
