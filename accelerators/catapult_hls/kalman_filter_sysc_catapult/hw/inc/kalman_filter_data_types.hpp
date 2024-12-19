// Copyright (c) 2011-2024 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0

#ifndef __DATATYPES__
#define __DATATYPES__

#include "ac_int.h"
#include "ac_fixed.h"
#include "kalman_filter_specs.hpp"
#include "ac_float.h"
#include <ac_std_float.h>
// #define FL_POINT

#define N const_mat_dim
#define FPDATA_WL DATA_WIDTH
#define FPDATA_WL DATA_WIDTH
#define FPDATA_IL DATA_WIDTH/2


typedef ac_int<DMA_WIDTH> DMA_WORD;
typedef ac_int<FPDATA_WL> FPDATA_WORD;
typedef ac_fixed<FPDATA_WL, FPDATA_IL> FPDATA;

#define __AC_FLOAT_ENABLE_ALPHA
// typedef ac_float<23, 1, 8> FLOAT_TYPE;

#ifdef FL_POINT
typedef float FLOAT_TYPE;
typedef FLOAT_TYPE FN_DATATYPE;
#else
typedef ac_ieee_float32 FLOAT_TYPE;
typedef FLOAT_TYPE FN_DATATYPE;
#endif

// inline void int2fp(const FPDATA_WORD& in, FLOAT_TYPE& out) {
//     float temp = *reinterpret_cast<const float*>(&in);  // reinterpret cast
//     out = FLOAT_TYPE(temp);  // Assign float to FLOAT_TYPE
// }

// inline void fp2int(const FLOAT_TYPE& in, FPDATA_WORD& out) {
//     float temp = in.to_float();  // Convert FLOAT_TYPE to float
//     out = *reinterpret_cast<FPDATA_WORD*>(&temp);  // reinterpret cast to FPDATA_WORD
// }

#ifdef FL_POINT
inline void int2fp(const FPDATA_WORD& in, FLOAT_TYPE& out) 
{ uint32_t data = in.to_uint(); float *ptr = (float *) &data; out = *ptr;\
}

inline void fp2int(const FLOAT_TYPE& in, FPDATA_WORD& out)
{ uint32_t *ptr = (uint32_t *) &in; out = *ptr; }

#else

// Function to convert FPDATA_WORD to FLOAT_TYPE
inline void int2fp(const FPDATA_WORD& in, FLOAT_TYPE& out) {
    out.set_data(in);
}

// Function to convert FLOAT_TYPE to FPDATA_WORD
inline void fp2int(const FLOAT_TYPE& in, FPDATA_WORD& out) {
    out = in.data(); // Fetch binary representation.
}

#endif
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
#ifdef FL_POINT
            // printf("%.30f ", matrix[i * cols + j]);
#else
            // printf("%.30f ", matrix[i * cols + j].to_float());
#endif
        }
        // cout << std::endl;
    }
}
#endif

inline void copymat(FN_DATATYPE A[N][N], FN_DATATYPE result[N][N], uint32_t kalman_mat_dim) {
    for (uint32_t i = 0; i < kalman_mat_dim; i++) {
        for (uint32_t j = 0; j < kalman_mat_dim; j++) {
            result[i][j] = A[i][j];
        }
    }
}

// inline void matrix_multiply(FN_DATATYPE* A, FN_DATATYPE* B, FN_DATATYPE* C, uint32_t n, uint32_t m, uint32_t p) {
inline void matrix_multiply(FN_DATATYPE A[MAX_MEAS_SIZE*MAX_MEAS_SIZE], FN_DATATYPE B[MAX_MEAS_SIZE*MAX_MEAS_SIZE], FN_DATATYPE C[MAX_MEAS_SIZE*MAX_MEAS_SIZE], uint32_t n, uint32_t m, uint32_t p) {
    // std::cout << "\tMatrix Multiply: " << n << "\t" << m << "\t" << p << "\n";;
    for (uint32_t i = 0; i < n; i++) {
        for (uint32_t j = 0; j < p; j++) {
            // C[i * p + j] = 0;
            C[i * p + j] = FN_DATATYPE(0.0);

            for (uint32_t k = 0; k < m; k++) {
                C[i * p + j] += A[i * m + k] * B[k * p + j];
            }
        }
    }
}



// inline void inverse_clean(FN_DATATYPE new_mat[TMP_MAX_SIZE_INV][TMP_MAX_SIZE_INV], FN_DATATYPE out[TMP_MAX_SIZE_INV][TMP_MAX_SIZE_INV], uint32_t meas_size_reg)
// {

// 		 FN_DATATYPE ratio;
// 		 uint32_t i,j,k;

// 		 if(meas_size_reg == 2){
// 			 FN_DATATYPE a = new_mat[0][0];
// 			 FN_DATATYPE b = new_mat[0][1];
// 			 FN_DATATYPE c = new_mat[1][0];
// 			 FN_DATATYPE d = new_mat[1][1];

// 			 FN_DATATYPE det = FN_DATATYPE((a * d) - (b * c));


// 			 if (det == FN_DATATYPE(0.0)) {
// 			     return;
// 			 }

// 			 out[0][0] = d / det;
// 			 out[0][1] = FN_DATATYPE(-1.0) * b / det;
// 			 out[1][0] = FN_DATATYPE(-1.0) * c / det;
// 			 out[1][1] = a / det;

// 			 return;
// 		 }

// 		 /* Applying Gauss Jordan Elimination */
// 		 for(i = 0; i < meas_size_reg; i++)
// 		 {
// 			  for(j = 0; j < meas_size_reg; j++)
// 			  {
// 				   if(i != j)
// 				   {
// 					    ratio = new_mat[j][i]/new_mat[i][i];
// 					    for(k = 0; k < meas_size_reg; k++)
// 					    {

// 					    	if(i == meas_size_reg-1){
// 					    		if(k == 0){//Calc the diagonal element first
// 					    			new_mat[j][j] = FN_DATATYPE(new_mat[j][j] - ratio*new_mat[i][j]);
// 					    		}
// 					    		else if(k == j){
// 					    			new_mat[j][0] = FN_DATATYPE(new_mat[j][0] - ratio*new_mat[i][0]);
// 					    		}
// 					    		else{
// 					    			new_mat[j][k] = FN_DATATYPE(new_mat[j][k] - ratio*new_mat[i][k]);
// 					    		}

// 					    		out[j][k] = FN_DATATYPE((out[j][k] - ratio*out[i][k]) / new_mat[j][j]);
// 					    	}
// 					    	else{

// 								new_mat[j][k] = FN_DATATYPE(new_mat[j][k] - ratio*new_mat[i][k]);

// 								if(i > 0)
// 									out[j][k] = FN_DATATYPE(out[j][k] - ratio*out[i][k]);
// 								else{ //(i == 0)
// 									if(i == k)
// 										out[i][k] = FN_DATATYPE(1);
// 									else
// 										out[i][k] = FN_DATATYPE(0);
// 									if(j == k){
// 										if(i == k)
// 											out[j][k] = FN_DATATYPE((FN_DATATYPE)1 - ratio);
// 											// out[j][k] = FN_DATATYPE(1 - ratio.to_float());
// 										else
// 											out[j][k] = FN_DATATYPE(1);
// 									}
// 									else{
// 										if(i == k)
// 											out[j][k] = -ratio;
// 										else
// 											out[j][k] = FN_DATATYPE(0);
// 									}
// 								}
// 					    	}

// 					    }

// 				   }

// 			  }
// 		 }
// 		 for(i = 0; i < meas_size_reg; i++)
// 		 {
// 			 out[meas_size_reg-1][i] = out[meas_size_reg-1][i] / new_mat[meas_size_reg-1][meas_size_reg-1];
// 		 }
// 		 return;
// }

// inline void inverse_clean(FN_DATATYPE* new_mat, FN_DATATYPE* out, uint32_t meas_size_reg)
inline void inverse_clean(FN_DATATYPE new_mat[TMP_MAX_SIZE_INV*TMP_MAX_SIZE_INV], FN_DATATYPE out[TMP_MAX_SIZE_INV*TMP_MAX_SIZE_INV], uint32_t meas_size_reg)
{
    FN_DATATYPE ratio;
    uint32_t i, j, k;

    if (meas_size_reg == 2) {
        FN_DATATYPE a = new_mat[0 * meas_size_reg + 0];
        FN_DATATYPE b = new_mat[0 * meas_size_reg + 1];
        FN_DATATYPE c = new_mat[1 * meas_size_reg + 0];
        FN_DATATYPE d = new_mat[1 * meas_size_reg + 1];

        FN_DATATYPE det = FN_DATATYPE((a * d) - (b * c));

        if (det == FN_DATATYPE(0.0)) {
            return;
        }

        out[0 * meas_size_reg + 0] = d / det;
        out[0 * meas_size_reg + 1] = FN_DATATYPE(-1.0) * b / det;
        out[1 * meas_size_reg + 0] = FN_DATATYPE(-1.0) * c / det;
        out[1 * meas_size_reg + 1] = a / det;

        return;
    }

    /* Applying Gauss Jordan Elimination */
    for (i = 0; i < meas_size_reg; i++) {
        for (j = 0; j < meas_size_reg; j++) {
            if (i != j) {
                ratio = new_mat[j * meas_size_reg + i] / new_mat[i * meas_size_reg + i];
                for (k = 0; k < meas_size_reg; k++) {
                    if (i == meas_size_reg - 1) {
                        if (k == 0) {
                            new_mat[j * meas_size_reg + j] -= ratio * new_mat[i * meas_size_reg + j];
                        } else if (k == j) {
                            new_mat[j * meas_size_reg + 0] -= ratio * new_mat[i * meas_size_reg + 0];
                        } else {
                            new_mat[j * meas_size_reg + k] -= ratio * new_mat[i * meas_size_reg + k];
                        }

                        out[j * meas_size_reg + k] =
                            (out[j * meas_size_reg + k] - ratio * out[i * meas_size_reg + k]) /
                            new_mat[j * meas_size_reg + j];
                    } else {
                        new_mat[j * meas_size_reg + k] -= ratio * new_mat[i * meas_size_reg + k];

                        if (i > 0) {
                            out[j * meas_size_reg + k] -= ratio * out[i * meas_size_reg + k];
                        } else {
                            if (i == k) {
                                out[i * meas_size_reg + k] = FN_DATATYPE(1);
                            } else {
                                out[i * meas_size_reg + k] = FN_DATATYPE(0);
                            }
                            if (j == k) {
                                if (i == k) {
                                    out[j * meas_size_reg + k] = FN_DATATYPE(1) - ratio;
                                } else {
                                    out[j * meas_size_reg + k] = FN_DATATYPE(1);
                                }
                            } else {
                                if (i == k) {
                                    out[j * meas_size_reg + k] = -ratio;
                                } else {
                                    out[j * meas_size_reg + k] = FN_DATATYPE(0);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    for (i = 0; i < meas_size_reg; i++) {
        out[(meas_size_reg - 1) * meas_size_reg + i] /=
            new_mat[(meas_size_reg - 1) * meas_size_reg + (meas_size_reg - 1)];
    }
    return;
}


// Utility function implementations

inline void gauss_inverse(FN_DATATYPE* A, FN_DATATYPE* A_inv, int n) {
    // Augmenting the matrix A with identity matrix of same dimensions
    FN_DATATYPE augmented[MAX_MEAS_SIZE * 2 * MAX_MEAS_SIZE];

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            augmented[i * 2 * n + j] = A[i * n + j];  // A portion
            augmented[i * 2 * n + (j + n)] = (i == j) ? FN_DATATYPE(1) : FN_DATATYPE(0);  // Identity portion
        }
    }

    for (int i = 0; i < n; i++) {
        int pivot_row = i;
        for (int j = i + 1; j < n; j++) {
            if (augmented[j * 2 * n + i] > augmented[pivot_row * 2 * n + i]) {
                pivot_row = j;
            }
        }

        if (pivot_row != i) {
            for (int k = 0; k < 2 * n; k++) {
                FN_DATATYPE temp = augmented[i * 2 * n + k];
                augmented[i * 2 * n + k] = augmented[pivot_row * 2 * n + k];
                augmented[pivot_row * 2 * n + k] = temp;
            }
        }

        FN_DATATYPE pivot = augmented[i * 2 * n + i];
        for (int k = 0; k < 2 * n; k++) {
            augmented[i * 2 * n + k] /= pivot;
        }

        for (int j = 0; j < n; j++) {
            if (j != i) {
                FN_DATATYPE factor = augmented[j * 2 * n + i];
                for (int k = 0; k < 2 * n; k++) {
                    augmented[j * 2 * n + k] -= factor * augmented[i * 2 * n + k];
                }
            }
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A_inv[i * n + j] = augmented[i * 2 * n + (j + n)];
        }
    }
}


// inline void matrix_add(FN_DATATYPE* A, FN_DATATYPE* B, FN_DATATYPE* C, int n, int m) {
inline void matrix_add(FN_DATATYPE A[MAX_MEAS_SIZE*MAX_MEAS_SIZE], FN_DATATYPE B[MAX_MEAS_SIZE*MAX_MEAS_SIZE], FN_DATATYPE C[MAX_MEAS_SIZE*MAX_MEAS_SIZE], int n, int m) {
    // std::cout << "\tMatrix Add: " << n << "\t" << m << "\n";;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            C[i * m + j] = FLOAT_TYPE(A[i * m + j] + B[i * m + j]);
        }
    }
}


// inline void matrix_transpose(FN_DATATYPE* A, FN_DATATYPE* AT, int n, int m) {
inline void matrix_transpose(FN_DATATYPE A[MAX_MEAS_SIZE*STATE_SIZE], FN_DATATYPE AT[MAX_MEAS_SIZE*STATE_SIZE], int n, int m) {
    // std::cout << "\tMatrix Transpose: " << n << "\t" << m << "\n";;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            AT[j * n + i] = A[i * m + j];
        }
    }
}

// inline void matrix_subtract(FN_DATATYPE* A, FN_DATATYPE* B, FN_DATATYPE* C, int n) {
inline void matrix_subtract(FN_DATATYPE A[MAX_MEAS_SIZE], FN_DATATYPE B[MAX_MEAS_SIZE], FN_DATATYPE C[MAX_MEAS_SIZE], int n) {
    // std::cout << "\tMatrix Subtract: " << n << "\n";;
    for (int i = 0; i < n; i++) {
        C[i] = FLOAT_TYPE(A[i] - B[i]);
    }
}

#endif
