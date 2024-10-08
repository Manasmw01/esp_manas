//Copyright (c) 2011-2024 Columbia University, System Level Design Group
//SPDX-License-Identifier: Apache-2.0

#include "testbench.hpp"
#include "ac_math/ac_random.h"
#include <mc_connections.h>
#include <mc_scverify.h>

#include <ac_float.h>

#include "A_array.h"
#include "H_array.h"
#include "initial_state_array.h"
#include "measurements_array.h"
#include "P_array.h"
#include "prediction_array.h"
#include "Q_array.h"
#include "real_array.h"
#include "W_array.h"


std::ofstream ofs;
std::ifstream ifs;

#define ERROR_THRESHOLD 0.1
int err=0;
FPDATA data;
float float_data;
FPDATA Pp_Final[STATE_SIZE][STATE_SIZE];

#ifdef GOLDEN_OP
inline void print_matrix_golden(float matrixx[N][N], uint32_t kalman_mat_rows)
{
    for (int i = 0; i < kalman_mat_rows; i++)
        for (int j = 0; j < kalman_mat_rows; j++)
            std::cout << std::setw(20) << matrixx[i][j] << ((j == kalman_mat_rows - 1) ? "\n": "\t");
    std::cout << "\n";
}

inline void print_vector_golden(float vec[N] , uint32_t kalman_mat_rows)
{
    for (int i = 0; i < kalman_mat_rows; i++)
            std::cout << std::setw(20) << vec[i] << "\t";
    std::cout << "\n";
}
inline void gauss_inverse_golden(float A[N][N], float inverse[N][N], uint32_t kalman_mat_rows) {
    // Augmenting the matrix A with identity matrix of same dimensions
    float augmented[N][2*N];
    for (int i = 0; i < kalman_mat_rows; i++) {
        for (int j = 0; j < kalman_mat_rows; j++) {
            augmented[i][j] = A[i][j];
            // Set the corresponding element in the augmented matrix to 1 if on the diagonal, 0 otherwise
            // augmented[i][j + N] = (i == j) ? 1.0 : 0.0;

            if(i == j)
                augmented[i][j + kalman_mat_rows] = 1;                
            else
                augmented[i][j + kalman_mat_rows] = 0;                
        }
    }

    // Applying Gauss-Jordan elimination
    for (int i = 0; i < kalman_mat_rows; i++) {
        // Find the pivot row and swap
        int pivot_row = i;
        for (int j = i + 1; j < kalman_mat_rows; j++) 
        {
            if (augmented[j][i] > augmented[pivot_row][i]) 
            {
                pivot_row = j;
            }
        }
        // Swap rows i and pivot_row
        if (pivot_row != i) 
        {
            for (int k = 0; k < 2*kalman_mat_rows; k++) 
            {
                float temp = augmented[i][k];
                // augmented[i][k] = augmented[pivot_row][k];
                float temp2 = augmented[pivot_row][k];
                augmented[i][k] = temp2;
                augmented[pivot_row][k] = temp;
            }
        }

        // Make the diagonal elements 1
        float pivot = augmented[i][i];
        for (int k = 0; k < 2*kalman_mat_rows; k++) 
        {
            float temp = augmented[i][k] / pivot;
            augmented[i][k] = temp;
        }

        // Make other elements in the column 0
        for (int j = 0; j < kalman_mat_rows; j++) 
        {
            if (j != i) 
            {
                float factor = augmented[j][i];
                for (int k = 0; k < 2*kalman_mat_rows; k++) {
                    // augmented[j][k] -= factor * augmented[i][k];
                    float temp = factor * augmented[i][k];
                    float temp2 = augmented[j][k] - temp;
                    augmented[j][k] = temp2;
                }
            }
        }
    }

    // Copy the inverse from the augmented matrix
    for (int i = 0; i < kalman_mat_rows; i++) 
    {
        for (int j = 0; j < kalman_mat_rows; j++) 
        {
            inverse[i][j] = augmented[i][j + kalman_mat_rows];
        }
    }
}


inline void multiplyMatrices_golden(float A[N][N], float B[N][N], float result[N][N], uint32_t kalman_mat_rows) {
    for (int i = 0; i < kalman_mat_rows; i++) {
        for (int j = 0; j < kalman_mat_rows; j++) {
            result[i][j] = 0;
            for (int k = 0; k < kalman_mat_rows; k++) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}


inline void multiplyMatrixVector_golden(float A[N][N], float vector[N], float result[N], uint32_t kalman_mat_rows) {
    for (int i = 0; i < kalman_mat_rows; i++) {
        result[i] = 0;
        for (int j = 0; j < kalman_mat_rows; j++) {
            result[i] += A[i][j] * vector[j];
        }
    }
}

inline void transposeMatrix_golden(float matrix[N][N], float result[N][N], uint32_t kalman_mat_rows) {
    for (int i = 0; i < kalman_mat_rows; i++) {
        for (int j = 0; j < kalman_mat_rows; j++) {
            result[j][i] = matrix[i][j];
        }
    }
}

inline void addMatrices_golden(float A[N][N], float B[N][N], float result[N][N], uint32_t kalman_mat_rows) {
    for (int i = 0; i < kalman_mat_rows; i++) {
        for (int j = 0; j < kalman_mat_rows; j++) {
            result[i][j] = A[i][j] + B[i][j];
        }
    }
}

inline void subtractMatrices_golden(float A[N][N], float B[N][N], float result[N][N], uint32_t kalman_mat_rows) {
    for (int i = 0; i < kalman_mat_rows; i++) {
        for (int j = 0; j < kalman_mat_rows; j++) {
            result[i][j] = A[i][j] - B[i][j];
        }
    }
}

inline void copymat_golden(float A[N][N], float result[N][N], uint32_t kalman_mat_rows) {
    for (int i = 0; i < kalman_mat_rows; i++) {
        for (int j = 0; j < kalman_mat_rows; j++) {
            result[i][j] = A[i][j];
        }
    }
}

inline void subtractVectors_golden(float vector1[N], float vector2[N], float result[N], uint32_t kalman_mat_rows) {
    for (int i = 0; i < kalman_mat_rows; i++) {
        result[i] = vector1[i] - vector2[i];
    }
}

inline void addVectors_golden(float vector1[N], float vector2[N], float result[N], uint32_t kalman_mat_rows) {
    for (int i = 0; i < kalman_mat_rows; i++) {
        result[i] = vector1[i] + vector2[i];
    }
}
#endif

void testbench::proc()
{
    conf_info.Reset();
    mac_n = 1;
    mac_vec = 100;
    mac_len = 64;
    kalman_iters = num_iterations; // Number of readings taken
    kalman_mat_rows = STATE_SIZE; // Number of num_iterations (matrix_dim: [X_GPS(i); X_pos(i); Y_GPS(i); Y_pos(i)])
    kalman_mat_cols = STATE_SIZE; // Number of num_iterations (matrix_dim: [X_GPS(i); X_pos(i); Y_GPS(i); Y_pos(i)])

    vec_X_address = 0;
    Mat_F_address = vec_X_address + (STATE_SIZE);
    Mat_Q_address = Mat_F_address + (STATE_SIZE * STATE_SIZE);
    Mat_R_address = Mat_Q_address + (STATE_SIZE * STATE_SIZE);
    Mat_H_address = Mat_R_address + (MEAS_SIZE * MEAS_SIZE);
    // vec_Z_address = Mat_H_address + (MEAS_SIZE * STATE_SIZE);
    Mat_P_address = Mat_H_address + (MEAS_SIZE * STATE_SIZE);
    constant_matrices_size = Mat_P_address + (STATE_SIZE * STATE_SIZE);

    measurement_vecs_base_address = constant_matrices_size;

    input_vecs_total_size = measurement_vecs_base_address + SAMPLES*MEAS_SIZE;

    output_size_per_iter = STATE_SIZE + STATE_SIZE*STATE_SIZE; // xp and Pp
    output_total_size = output_size_per_iter;

    wait();

    CCS_LOG("--------------------------------");
    CCS_LOG("ESP - kalman_filter [Catapult HLS SystemC]");
    CCS_LOG("--------------------------------");

#if (DMA_WORD_PER_BEAT == 0)
    in_words_adj = input_vecs_total_size;
    out_words_adj = output_total_size;
#else
    in_words_adj = round_up(input_vecs_total_size, DMA_WORD_PER_BEAT);
    out_words_adj = round_up(output_total_size, DMA_WORD_PER_BEAT);
#endif
    in_size = in_words_adj * (mac_n);
    out_size = out_words_adj * (mac_n);

    CCS_LOG( "  - DMA width: "<< DMA_WIDTH);
    CCS_LOG( "  - DATA width: "<< DATA_WIDTH);
    CCS_LOG("--------------------------------");


    in = new ac_int<DATA_WIDTH,false>[in_size];
    in_float = new float[in_size];
    gold= new ac_int<DATA_WIDTH,false>[out_size];
    gold_float = new float[out_size];
    master_array = new float[in_size];
    golden_array = new float[out_size];


    std::cout << "vec_X_address\t" << vec_X_address << std::endl;
    std::cout << "Mat_F_address\t" << Mat_F_address << std::endl;
    std::cout << "Mat_Q_address\t" << Mat_Q_address << std::endl;
    std::cout << "Mat_R_address\t" << Mat_R_address << std::endl;
    std::cout << "Mat_H_address\t" << Mat_H_address << std::endl;
    std::cout << "Mat_P_address\t" << Mat_P_address << std::endl;

    std::cout << "constant_matrices_size\t" << constant_matrices_size << std::endl << std::endl;
    std::cout << "input_vecs_total_size\t" << input_vecs_total_size << std::endl;
    std::cout << "output_size_per_iter\t" << output_size_per_iter << std::endl;
    std::cout << "output_total_size\t" << output_total_size << std::endl;
    std::cout << "in_size\t" << in_size << std::endl;
    std::cout << "out_size\t" << out_size << std::endl;

    partition();
    std::cout << "Partition completed\n";

    // print_variables();
    // std::cout << "Print variables completed\n";

    single_input_array(); // Updates in_float by merging all the inputs
    std::cout << "Single input array completed\n";

    std::cout << "load_data\t" << input_vecs_total_size << std::endl;
    load_data(in_float, input_vecs_total_size); // Writes data in mem[i]
    std::cout << "Load datafloat done\n";

    do_config();
    std::cout << "Do config done\n";

    dump_memory();
    std::cout << "Dump memory completed\n";

    validate();

    
    sc_stop();

    wait();
}


void testbench::load_data(float *inn, uint32_t inn_size)
{
    for (uint32_t i = 0; i < inn_size / DMA_WORD_PER_BEAT; i++)  {
        ac_int<DMA_WIDTH> data_bv;
        for (int wordd = 0; wordd < DMA_WORD_PER_BEAT; wordd++)
        {
            // ac_ieee_float32 data = inn[i* DMA_WORD_PER_BEAT + wordd];
            // ac_float< 5, 3, 3, AC_RND> data = inn[i* DMA_WORD_PER_BEAT + wordd];
            // cout << "OUT 1" << endl; 
            ac_float< DATA_WIDTH, FPDATA_IL, 8, AC_RND> data = inn[i* DMA_WORD_PER_BEAT + wordd];
            // cout << "OUT 2" << endl; 

            // FPDATA fpdata=data.convert_to_ac_fixed<FPDATA_WL,FPDATA_IL,true,AC_TRN, AC_WRAP>();
            // FPDATA fpdata=static_cast<FPDATA>(data.to_double());
            FPDATA fpdata=data.to_ac_fixed();
            // if(i < 20)
            // cout << "ac_float[" <<  i << "]:" << std::setprecision(20) << fpdata << endl; 

            
            FPDATA_WORD fpdata_word;
            fpdata_word.set_slc(0,fpdata.slc<FPDATA_WL>(0));

            data_bv.set_slc(wordd*FPDATA_WL,fpdata_word);
        }
        mem[i] = data_bv;
    }
}

void testbench::partition()
{
    for (int i = 0; i < STATE_SIZE; i++) 
    {
        master_array[vec_X_address + i] = initial[i];
        // cout << "vecx[" << vec_X_address << "]:" << initial[i] << endl; 
    }

    for (int i = 0; i < STATE_SIZE; i++) {
        for (int j = 0; j < STATE_SIZE; j++) {
            master_array[Mat_F_address + i * STATE_SIZE + j] = A[i * STATE_SIZE + j];
            // cout << "Mat_F[" <<  Mat_F_address + i * STATE_SIZE + j << "]:" << std::setprecision(20) << A[i * STATE_SIZE + j] << endl; 
        }
    }

    for (int i = 0; i < STATE_SIZE; i++) {
        for (int j = 0; j < STATE_SIZE; j++) {
            master_array[Mat_Q_address + i * STATE_SIZE + j] = W[i * STATE_SIZE + j];
        }
    }

    for (int i = 0; i < MEAS_SIZE; i++) {
        for (int j = 0; j < MEAS_SIZE; j++) {
            master_array[Mat_R_address + i * MEAS_SIZE + j] = Q[i * MEAS_SIZE + j];
        }
    }

    for (int i = 0; i < STATE_SIZE; i++) {
        for (int j = 0; j < MEAS_SIZE; j++) {
            master_array[Mat_H_address + i * MEAS_SIZE + j] = H[i * MEAS_SIZE + j];
        }
    }

    for (int i = 0; i < STATE_SIZE; i++)
    {
        for (int j = 0; j < STATE_SIZE; j++)
        {
            master_array[Mat_P_address + (i*STATE_SIZE + j)] = 0; 
            // std::cout << "Mat_P[" << Mat_P_address + (i*STATE_SIZE + j) << "]:\t" << master_array[Mat_P_address + (i*STATE_SIZE + j)] << std::endl;
        }
    }
    for(int iter = 1; iter <= SAMPLES; iter++)
    {
        for(int i = MEAS_SIZE*iter; i < MEAS_SIZE*(iter+1); i++)
        {
            master_array[measurement_vecs_base_address + i-MEAS_SIZE*iter + MEAS_SIZE*(iter-1)] = measurements[i];    
            // if(i == MEAS_SIZE*(iter+1) - 1)
            // if(i < (MEAS_SIZE*iter + 6))
            // {
            //     std::cout << "iter: " << iter << "\ti:" << i << "\tMEAS_SIZE*iter: " << MEAS_SIZE*iter << "\tTB ZK[" << measurement_vecs_base_address + i-MEAS_SIZE*iter + MEAS_SIZE*(iter-1) << "]:" << master_array[measurement_vecs_base_address + i-MEAS_SIZE*iter + MEAS_SIZE*(iter-1)] << "\t" << measurements[i] << std::endl;
            // }
                // std::cout << "TB ZK:\t" << measurement_vecs_base_address + i-MEAS_SIZE*iter << std::endl;

        }
        std::cout << std::endl;
    }
}





void testbench::do_config()
{
    {
        conf_info_t config;

        /* <<--params-->> */
        config.mac_n = mac_n;
        config.mac_vec = mac_vec;
        config.mac_len = mac_len;
        config.kalman_iters = kalman_iters;
        config.kalman_mat_rows = kalman_mat_rows;

        config.constant_matrices_size = constant_matrices_size;
        config.measurement_vecs_base_address = measurement_vecs_base_address;

        config.vec_X_address = vec_X_address;
        config.Mat_F_address = Mat_F_address;
        config.Mat_Q_address = Mat_Q_address;
        config.Mat_R_address = Mat_R_address;
        config.Mat_H_address = Mat_H_address;
        // config.vec_Z_address = vec_Z_address;
        config.Mat_P_address = Mat_P_address;

        config.input_vecs_total_size = input_vecs_total_size;
        config.output_total_size = output_total_size;

        conf_info.Push(config);
    }
}

void testbench::dump_memory()
{
    std::cout << "Entered Dump Memory\n";
    do 
    {
        wait(); 
    } while (!acc_done.read());
    int offset=in_size;
    std::cout << "Dump memory offset address: "<< offset << "\n";

    out= new ac_int<DATA_WIDTH,false>[out_size];
    out_float= new float[out_size];

    offset = offset / DMA_WORD_PER_BEAT;
    ofs.open("accelerator_output.txt", std::ofstream::out);
    for (uint32_t i = 0; i < out_size / DMA_WORD_PER_BEAT; i++)
    {
        for (uint32_t wordd = 0; wordd < DMA_WORD_PER_BEAT; wordd++)
        {
            out[i * DMA_WORD_PER_BEAT + wordd] = mem[offset + i].slc<DATA_WIDTH>(wordd*DATA_WIDTH);
            FPDATA out_fixed = 0;
            int2fx(out[i * DMA_WORD_PER_BEAT + wordd],out_fixed);
            if(i >= out_size - (STATE_SIZE*STATE_SIZE))
            {
                ofs << out_fixed << std::endl;
            }

            ofs << i << ": " << out[i * DMA_WORD_PER_BEAT + wordd] << std::endl;
        }
    }
    ofs.close();
}


void testbench::validate()
{
    int tot_errors = 0;

    for (uint32_t i = 0; i < mac_n; i++)
        for (uint32_t j = 0; j < output_size_per_iter*num_iterations; j++)
        {

            // FPDATA out_gold_fx = 0;
            // int2fx(gold[i * out_words_adj + j],out_gold_fx);

            FPDATA out_res_fx = 0;
            int2fx(out[i * out_words_adj + j],out_res_fx);

// cout << "\nTESTTT\n";
// cout << output_size_per_iter*num_iterations << "\n";
// cout << output_size_per_iter*num_iterations - output_size_per_iter << "\n";
            if(j >= output_size_per_iter*num_iterations - output_size_per_iter)
            {
                    if(j%output_size_per_iter == 0)
                        std::cout << "\nXp_design[" << j/output_size_per_iter << "]:\t";
                    else if (j%output_size_per_iter == STATE_SIZE)
                    {
                        std::cout << "\nPp_design[" << j/output_size_per_iter << "]:\t";                        
                    }
                    if (j%STATE_SIZE == 0)
                    {
                        std::cout << "\n";                        
                    }
                std::cout << std::setprecision(20)  << out_res_fx << "\t";
            }
        }
    // float accelerator_validate_array[100];
    // float golden_validate_array[100];

    // FILE *file;
    // char line[100]; // Adjust the size as per your needs
    // // Open the file in read mode
    // file = fopen("golden_output.txt", "r");
    // if (file == NULL) {
    //     printf("Error opening the file.\n");
    // }

    // // Read each line from the file
    // uint32_t ind = 0;
    // while (fgets(line, sizeof(line), file)) {
    //     // Convert the line to a float
    //     float value = strtof(line, NULL);
    //     // Check if the conversion was successful
    //     if (value != 0.0f || (value == 0.0f && line[0] == '0')) {
    //         golden_validate_array[ind] = value;
    //         // Print the float value
    //         // printf("%f\n", value);
    //     } 
    //     ind++;
    // }
    // // Close the file
    // fclose(file);


    // // Open the file in read mode
    // file = fopen("accelerator_output.txt", "r");
    // if (file == NULL) {
    //     printf("Error opening the file.\n");
    // }

    // // Read each line from the file
    // ind = 0;
    // while (fgets(line, sizeof(line), file)) {
    //     // Convert the line to a float
    //     float value = strtof(line, NULL);
    //     // Check if the conversion was successful
    //     if (value != 0.0f || (value == 0.0f && line[0] == '0')) {
    //         accelerator_validate_array[ind] = value;
    //         // Print the float value
    //         // printf("%f\n", value);
    //     } 
    //     ind++;
    // }
    // // Close the file
    // fclose(file);

    // for(int i = 0; i < STATE_SIZE*STATE_SIZE; i++)
    // {
    //     // cout << accelerator_validate_array[i] << "\t" << golden_validate_array[i] << "\n";
    //       if (accelerator_validate_array[i] != golden_validate_array[i])
    //     {
    //         float MSE = (accelerator_validate_array[i]-golden_validate_array[i])*(accelerator_validate_array[i]-golden_validate_array[i])
    //             / golden_validate_array[i];

    //         if (MSE > ERROR_THRESHOLD)
    //         {
    //             printf("L2: output[%d] = %f (expected: %f)",
    //                         i, accelerator_validate_array[i], golden_validate_array[i]);
    //             tot_errors += 1;
    //         }

       
    //     }
    // }
    // if (tot_errors == 0)
    // {
    //     printf("------------------------------------\n");
    //     printf("  Validation succeeded!  \n");
    //     printf("------------------------------------\n");
    // }
    // else
    // {
    //     printf("------------------------------------\n");
    //     printf("  Validation failed!  \n");
    //     printf("------------------------------------\n");
    // }     
}

void testbench::single_input_array()
{
    for (uint32_t i= 0; i< mac_n ; i++)
    {
        for (uint32_t j=0; j< in_size ; j+=1)
        {
            float_data = master_array[j];
            FPDATA_WORD data_int32;
            data_int32.set_slc(0,data.slc<DATA_WIDTH>(0));
            in_float[i*in_words_adj+j]=float_data;
        }
    }
}