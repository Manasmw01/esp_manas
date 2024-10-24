// Copyright (c) 2011-2023 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0

#ifndef __FUNCTIONS_HPP__
#define __FUNCTIONS_HPP__

#include "kalman_filter.hpp"
void kalman_filter_sysc_catapult::compute_req(uint32_t iter, uint32_t kalman_iters, uint32_t kalman_mat_dim, uint32_t constant_matrices_size, bool pingpong, bool out_pingpong)
{
    
        for (uint32_t i = 0; i < MEAS_SIZE; i++)
        {
            plm_RRq<in_as,inrp> rreq;
            rreq.indx[0]=i;
            if(pingpong)            
                in_ping_ra.Push(rreq);
            else
                in_pong_ra.Push(rreq);
        }
        
        for (uint32_t i = 0; i < constant_matrices_size; i++)
        {
            plm_RRq<in_as,inrp> rreq;
            rreq.indx[0]=i;
            if(pingpong)            
                in_b_ping_ra.Push(rreq);
            else
                in_b_pong_ra.Push(rreq);
        }

        if(iter > 0)
        {
            uint32_t previous_state_elements = kalman_mat_dim + (kalman_mat_dim * kalman_mat_dim);    
            for (uint32_t k = 0; k < previous_state_elements; k++)
            {
                plm_RRq<out_as,outrp> rreq;
                rreq.indx[0]=k;
                if (pingpong)
                    xp_ping_ra.Push(rreq);
                else
                    xp_pong_ra.Push(rreq);
            }
        }
        wait();
        //keep
        sync_comp.sync_out();
}


void kalman_filter_sysc_catapult::compute(uint32_t iter, uint32_t kalman_iters, uint32_t kalman_mat_dim, 
                                    uint32_t vec_X_address, uint32_t Mat_F_address  , uint32_t Mat_Q_address, uint32_t Mat_R_address, 
                                    uint32_t Mat_H_address, uint32_t Mat_P_address, uint32_t constant_matrices_size, bool pingpong, bool out_pingpong)
{

    FLOAT_TYPE vec_X[STATE_SIZE];
    FLOAT_TYPE Mat_P[STATE_SIZE * STATE_SIZE];
    FLOAT_TYPE Mat_F[STATE_SIZE * STATE_SIZE];
    FLOAT_TYPE Mat_Q[STATE_SIZE * STATE_SIZE];
    FLOAT_TYPE Mat_R[MEAS_SIZE * MEAS_SIZE];
    FLOAT_TYPE Mat_H[MEAS_SIZE * STATE_SIZE];
    FLOAT_TYPE vec_Z[MEAS_SIZE];

    FLOAT_TYPE Mat_K[STATE_SIZE * MEAS_SIZE];
    FLOAT_TYPE Mat_I[STATE_SIZE*STATE_SIZE];
    
    
    FLOAT_TYPE vec_Z_fp[MEAS_SIZE];  // Floating-point version of vec_Z
    FLOAT_TYPE  Pp_cpp[const_mat_dim][const_mat_dim];
    FLOAT_TYPE  X_cpp[const_mat_dim];
    FLOAT_TYPE  output_to_send[const_mat_dim];

    FPDATA_WORD important_matrices_word;
    FPDATA input_regs_fx;
    FLOAT_TYPE input_regs_fp;

    // FPDATA_WORD input_measurements_word[200];
    // FPDATA input_measurements_fx[200];

    FPDATA_WORD input_measurements_word;
    FLOAT_TYPE input_measurements_fx;

    FLOAT_TYPE zk[const_mat_dim];

    // uint32_t current_address = 0;


    for (uint32_t k = 0; k < MEAS_SIZE; k++)
    {
        if(pingpong)            
            input_measurements_word=in_ping_rd.Pop().data[0];
        else
            input_measurements_word=in_pong_rd.Pop().data[0];            

        // int2fx(input_measurements_word,input_measurements_fx);
            // vec_Z[k] = input_measurements_fx;
        FLOAT_TYPE recovered_float;
        int2fp(input_measurements_word, recovered_float);
            vec_Z[k] = recovered_float;
        // std::cout << "Vec_Z FLOAT_TYPE: " << recovered_float.to_double() << std::endl;

    }

    uint32_t rows = 0;
    uint32_t cols = 0;
    // cout << "functions compute:(" << constant_matrices_size<< "):\t" << endl; 

    for (uint32_t current_address = 0; current_address < constant_matrices_size;  current_address++)
    {
        if(pingpong)            
            important_matrices_word=in_b_ping_rd.Pop().data[0];
        else
            important_matrices_word=in_b_pong_rd.Pop().data[0];

        // int2fx(important_matrices_word,input_regs_fx);
        int2fp(important_matrices_word,input_regs_fp);
        if(current_address >= vec_X_address && current_address < Mat_F_address)
        {
            output_to_send[current_address - vec_X_address] = input_regs_fp;
            vec_X[current_address - vec_X_address] = input_regs_fp;
            // cout << "vecX[" <<  current_address << "]\t" << std::setprecision(20) << vec_X[current_address - vec_X_address] << "\n";           
        }
        if(current_address >= Mat_F_address && current_address < Mat_Q_address)
        {
            Mat_F[current_address - Mat_F_address] = input_regs_fp;
            // cout << "Mat_F[" << current_address << "]\t" << std::setprecision(20) << Mat_F[current_address - Mat_F_address] << "\n";           
        }

        if(current_address >= Mat_Q_address && current_address < Mat_R_address)
        {
            Mat_Q[current_address - Mat_Q_address] = input_regs_fp;
            // cout << "Mat_Q[" << current_address << "]\t" << std::setprecision(20) << Mat_Q[current_address - Mat_Q_address] << "\n";           
        }

        if(current_address >= Mat_R_address && current_address < Mat_H_address)
        {
            Mat_R[current_address - Mat_R_address] = input_regs_fp;
            // cout << "Mat_R[" << current_address << "]\t" << std::setprecision(20) << Mat_R[current_address - Mat_R_address] << "\n";
        }

        if(current_address >= Mat_H_address && current_address < Mat_P_address)
        {
            Mat_H[current_address - Mat_H_address] = input_regs_fp;
            // cout << "Mat_H[" << current_address << "]\t" << std::setprecision(20) << Mat_H[current_address - Mat_H_address] << "\n";
        }
        if(iter == 0)
            if(current_address >= Mat_P_address && current_address < constant_matrices_size)
            {
                Mat_P[current_address - Mat_P_address] = input_regs_fp;
                // cout << "Mat_P[" << current_address << "]\t" << std::setprecision(20) << Mat_P[current_address - Mat_P_address] << "\n";
            }   
    }


    FPDATA_WORD op_word;
    FPDATA op_fx;
    FLOAT_TYPE op_fp;
    
    if(iter > 0)
    {
        int indexxxx = 0;
        int rows = 0;
        int cols = 0;
        uint32_t previous_states_elements = kalman_mat_dim + (kalman_mat_dim * kalman_mat_dim);    
        // for (uint32_t k = 0; k < kalman_mat_dim; k++)
        for (uint32_t k = 0; k < previous_states_elements; k++)
        {
            if (pingpong)
                op_word=xp_ping_rd.Pop().data[0];
            else
                op_word=xp_pong_rd.Pop().data[0];
            // int2fx(op_word,op_fx);
            int2fp(op_word,op_fp);

            // floating op_float;
            // FPDATA op_fx_2;
            // op_fx_2.set_slc(0, op_word.template slc<FPDATA_WL>(0));
            // op_float = floating(op_fx_2);


            if(k < kalman_mat_dim)
            {
                vec_X[k] = op_fp; 
                // std::cout << "vec_X" << std::endl;
                // std::cout << "\t" << std::setprecision(20) << vec_X[k] << std::endl;
                // std::cout << "\t" << std::setprecision(20) << op_float << std::endl;
            }

            if(k >= kalman_mat_dim)
            {
                Mat_P[indexxxx] = op_fp;
                indexxxx++;
            }
        }

    // printf("(%d), NEW Mat_P\n", iter);
    // print_matrix_new(Mat_P, STATE_SIZE, STATE_SIZE);

    // printf("(%d), NEW vec_X\n", iter);
    // print_matrix_new(vec_X, STATE_SIZE, 1);

    }

    for (int i = 0; i < STATE_SIZE; i++) 
    {
        for (int j = 0; j < STATE_SIZE; j++) 
        {
            Mat_I[i * STATE_SIZE + j] = (i == j) ? 1.0 : 0.0;
        }
    }
    // printf("(%d), START vec_Z\n", iter);
    // print_matrix_new(vec_Z, 1, MEAS_SIZE);

    // printf("(%d), START Mat_P\n", iter);
    // print_matrix_new(Mat_P, STATE_SIZE, STATE_SIZE);

    // printf("(%d), START vec_X\n", iter);
    // print_matrix_new(vec_X, STATE_SIZE, 1);
        // #ifdef PRINT_STATEMENTS
        // std::cout << "PHI: " << iter << "\n";
        // print_matrix(phi_cpp, kalman_mat_dim);
        // #endif

        // #ifdef PRINT_STATEMENTS
        // std::cout << "X_Cpp: " << iter << "\n";
        // print_vector(X_cpp, kalman_mat_dim);
    // #endif
        FLOAT_TYPE Y[MEAS_SIZE];
        FLOAT_TYPE Mat_S[MEAS_SIZE*MEAS_SIZE];
        FLOAT_TYPE Mat_S_2D[MEAS_SIZE][MEAS_SIZE];


    // FPDATA HtF[MEAS_SIZE * STATE_SIZE];
    FLOAT_TYPE HtF[MEAS_SIZE * MEAS_SIZE];
    matrix_multiply(Mat_H, Mat_F, HtF, MEAS_SIZE,   STATE_SIZE, STATE_SIZE); // xp = A*x2

    // printf("HtF\n");
    // print_matrix_new(HtF, MEAS_SIZE, STATE_SIZE);

    FLOAT_TYPE H_F_X[MEAS_SIZE];
    matrix_multiply(HtF, vec_X, H_F_X, MEAS_SIZE,   STATE_SIZE, 1); // xp = A*x2
    // printf("H_F_X\n");
    // print_matrix_new(H_F_X, MEAS_SIZE, 1);

    // std::cout << "LINE 231" << std::endl;

    matrix_subtract(vec_Z, H_F_X, Y, MEAS_SIZE); // P3 = Pp - (K * H * Pp)
    // printf("Y\n");
    // print_matrix_new(Y, MEAS_SIZE, 1);

    FLOAT_TYPE F_Transpose[STATE_SIZE * STATE_SIZE];
    matrix_transpose(Mat_F, F_Transpose, STATE_SIZE, STATE_SIZE); // A^T
    // printf("F_Transpose\n");
    // print_matrix_new(F_Transpose, STATE_SIZE, STATE_SIZE);

    FLOAT_TYPE H_Transpose[STATE_SIZE * MEAS_SIZE];
    matrix_transpose(Mat_H, H_Transpose, MEAS_SIZE, STATE_SIZE); // A^T
    // printf("H_Transpose\n");
    // print_matrix_new(H_Transpose, STATE_SIZE, MEAS_SIZE);

    FLOAT_TYPE FtP[STATE_SIZE * STATE_SIZE];
    matrix_multiply(Mat_F, Mat_P, FtP, STATE_SIZE,   STATE_SIZE, STATE_SIZE); 
    // printf("FtP\n");
    // print_matrix_new(FtP, STATE_SIZE, STATE_SIZE);

    FLOAT_TYPE F_P_FT[STATE_SIZE * STATE_SIZE];
    matrix_multiply(FtP, F_Transpose, F_P_FT, STATE_SIZE,   STATE_SIZE, STATE_SIZE); 
    // printf("F_P_FT\n");
    // print_matrix_new(F_P_FT, STATE_SIZE, STATE_SIZE);

    FLOAT_TYPE F_P_FT_Q[STATE_SIZE * STATE_SIZE];
    matrix_add(F_P_FT, Mat_Q, F_P_FT_Q, STATE_SIZE, STATE_SIZE); // Pp = (A * P2) * A^T + Q
    // printf("F_P_FT_Q\n");
    // print_matrix_new(F_P_FT_Q, STATE_SIZE, STATE_SIZE);


    FLOAT_TYPE F_P_FT_Q_times_HT[STATE_SIZE*MEAS_SIZE];
    matrix_multiply(F_P_FT_Q, H_Transpose, F_P_FT_Q_times_HT, STATE_SIZE,   STATE_SIZE, MEAS_SIZE); 
    // printf("F_P_FT_Q_times_HT\n");
    // print_matrix_new(F_P_FT_Q_times_HT, STATE_SIZE, MEAS_SIZE);

    FLOAT_TYPE H_times_F_P_FT_Q_times_HT[MEAS_SIZE * MEAS_SIZE];
    matrix_multiply(Mat_H, F_P_FT_Q_times_HT, H_times_F_P_FT_Q_times_HT, MEAS_SIZE,   STATE_SIZE, MEAS_SIZE); 
    // printf("H_times_F_P_FT_Q_times_HT\n");
    // print_matrix_new(H_times_F_P_FT_Q_times_HT, MEAS_SIZE, MEAS_SIZE);
    matrix_add(H_times_F_P_FT_Q_times_HT, Mat_R, Mat_S, MEAS_SIZE, MEAS_SIZE); // Pp = (A * P2) * A^T + Q
    // printf("Mat_S\n");
    // print_matrix_new(Mat_S, MEAS_SIZE, MEAS_SIZE);

    FLOAT_TYPE S_inv[MEAS_SIZE * MEAS_SIZE];
    FLOAT_TYPE S_inv_2D[MEAS_SIZE][MEAS_SIZE];
    // gauss_inverse(Mat_S, S_inv, MEAS_SIZE); 

    // // INVERSE CLEAN IMPLEMENTATION STARTS
    for (int i = 0; i < MEAS_SIZE; i++)  
    {
        for (int j = 0; j < MEAS_SIZE; j++) 
        {
            Mat_S_2D[i][j] = Mat_S[i * MEAS_SIZE + j]; // Accessing the 1D array using the row-major order formula
        }
    }
    inverse_clean(Mat_S_2D, S_inv_2D);
    for (int i = 0; i < MEAS_SIZE; i++)  
    {
        for (int j = 0; j < MEAS_SIZE; j++) 
        {
            S_inv[i * MEAS_SIZE + j] = S_inv_2D[i][j];  // Converting 2D element back to 1D
        }
    }
    // printf("S_inv\n");
    // print_matrix_new(S_inv, MEAS_SIZE, MEAS_SIZE);

    matrix_multiply(F_P_FT_Q_times_HT, S_inv, Mat_K, STATE_SIZE,   MEAS_SIZE, MEAS_SIZE); 
    // printf("Mat_K\n");
    // print_matrix_new(Mat_K, STATE_SIZE, MEAS_SIZE);


    FLOAT_TYPE FtX[STATE_SIZE];
    matrix_multiply(Mat_F, vec_X, FtX, STATE_SIZE,   STATE_SIZE, 1); 
    // printf("FtX\n");
    // print_matrix_new(FtX, STATE_SIZE, 1);

    // FPDATA Mat_K[STATE_SIZE * MEAS_SIZE];
    // FPDATA Y[MEAS_SIZE];
    FLOAT_TYPE KtH[STATE_SIZE * MEAS_SIZE]; // Added Oct 8

    FLOAT_TYPE KtY[STATE_SIZE];
    matrix_multiply(Mat_K, Y, KtY, STATE_SIZE,   MEAS_SIZE, 1); 
    // printf("KtY\n");
    // print_matrix_new(KtY, STATE_SIZE, 1);

    matrix_add(FtX, KtY, vec_X, STATE_SIZE, 1); // Pp = (A * P2) * A^T + Q
    // printf("vec_X\n");
    // print_matrix_new(KtY, STATE_SIZE, 1);

    matrix_multiply(Mat_K, Mat_H, KtH, STATE_SIZE,   MEAS_SIZE, STATE_SIZE); 
    // printf("KtH\n");
    // print_matrix_new(KtH, STATE_SIZE, STATE_SIZE);
    // std::cout << "WORKING TILL HERE" << std::endl;

    FLOAT_TYPE I_minus_KtH[STATE_SIZE * STATE_SIZE];
    matrix_subtract(Mat_I, KtH, I_minus_KtH, STATE_SIZE*STATE_SIZE); // P3 = Pp - (K * H * Pp)

    // printf("I_minus_KtH\n");
    // print_matrix_new(I_minus_KtH, STATE_SIZE, STATE_SIZE);

    matrix_multiply(I_minus_KtH, F_P_FT_Q, Mat_P, STATE_SIZE,   STATE_SIZE, STATE_SIZE);
    // printf("Mat_P\n");
    // print_matrix_new(Mat_P, STATE_SIZE, STATE_SIZE);

    FLOAT_TYPE send_output_array[STATE_SIZE + (STATE_SIZE * STATE_SIZE)];
    
    uint32_t input_ptr;
    input_ptr = 0;
    uint32_t output_elements = kalman_mat_dim + (kalman_mat_dim * kalman_mat_dim);

    // printf("(%d), Updating Mat_P\n", iter);
    // print_matrix_new(Mat_P, STATE_SIZE, STATE_SIZE);

    // printf("(%d), Updating vec_X\n", iter);
    // print_matrix_new(vec_X, STATE_SIZE, 1);

    // std::cout << "(" << iter << "): Predicted vec_X\t";
    // for (int i = 0; i < STATE_SIZE; i++) 
    // {
    //     std::cout << std::setprecision(20) << vec_X[i] << "\t";
    // }
    // std::cout << "\n";

        for (uint32_t k = 0; k < STATE_SIZE; k++)
            send_output_array[k] = vec_X[k];
            // send_output_array[k] = output_to_send[k];
            // send_output_array[k] = Xp[k];

        for (int i = 0; i < STATE_SIZE; i++) {
            for (int j = 0; j < STATE_SIZE; j++) {
                send_output_array[STATE_SIZE + i * STATE_SIZE + j] = Mat_P[i * STATE_SIZE + j];
            }
        }

        // for (uint32_t i = 0; i < kalman_mat_dim; i++)
        //     for (uint32_t j = 0; j < kalman_mat_dim; j++)
        //         send_output_array[kalman_mat_dim + (i*kalman_mat_dim + j)] = Pp_cpp[i][j]; 

        for (uint32_t k = 0; k < output_elements; k++)
        {
            FPDATA_WORD out_word;
            fp2int(send_output_array[k], out_word);

            plm_WR<out_as,outwp> wreq;
            wreq.data[0]=out_word;
            wreq.indx[0]=input_ptr;

            if (out_pingpong)
                out_ping_w.Push(wreq);
            else
                out_pong_w.Push(wreq);
            input_ptr++;
        }
        input_ptr = 0;
        FLOAT_TYPE next_state_arr_fx[STATE_SIZE + (STATE_SIZE * STATE_SIZE)];
        for (uint32_t k = 0; k < kalman_mat_dim; k++)
            next_state_arr_fx[k] = vec_X[k];
            // next_state_arr_fx[k] = X_cpp[k];
            

        for (int i = 0; i < STATE_SIZE; i++) {
            for (int j = 0; j < STATE_SIZE; j++) {
                next_state_arr_fx[STATE_SIZE + i * STATE_SIZE + j] = Mat_P[i * STATE_SIZE + j];
            }
        }
        // for (uint32_t i = 0; i < kalman_mat_dim; i++)
        //     for (uint32_t j = 0; j < kalman_mat_dim; j++)
        //         next_state_arr_fx[kalman_mat_dim + (i*kalman_mat_dim + j)] = Pp_cpp[i][j]; 

        uint32_t next_state_size = kalman_mat_dim + (kalman_mat_dim * kalman_mat_dim);
        for (uint32_t k = 0; k < next_state_size; k++)
        {
            // std::cout << "k:" << k << "\n";
            FPDATA_WORD xp_word;
            fp2int(next_state_arr_fx[k], xp_word);

            plm_WR<out_as,outwp> wreq1;
            wreq1.data[0]=xp_word;
            wreq1.indx[0]=input_ptr;

            if (out_pingpong)
                xp_ping_w.Push(wreq1);
            else
                xp_pong_w.Push(wreq1);
            input_ptr++;
        }
    wait();
    // //keep
    sync_comp.sync_in();
}

void plusone(FPDATA& regg) 
{
    regg = regg + 1;
}
#endif // __FUNCTIONS_HPP__