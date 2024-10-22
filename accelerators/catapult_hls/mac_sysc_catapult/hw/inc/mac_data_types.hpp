// Copyright (c) 2011-2024 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0

#ifndef __DATATYPES__
#define __DATATYPES__

#include "ac_int.h"
#include "ac_fixed.h"
#include "mac_specs.hpp"

#define FPDATA_WL DATA_WIDTH
#define FPDATA_IL DATA_WIDTH/2


typedef ac_int<DMA_WIDTH> DMA_WORD;
typedef ac_int<FPDATA_WL> FPDATA_WORD;
typedef ac_fixed<FPDATA_WL, FPDATA_IL> FPDATA;


typedef ac_float<23, 0, 8> FLOAT_TYPE;


// Function to convert FPDATA_WORD to FLOAT_TYPE
inline void int2fp(const FPDATA_WORD& in, FLOAT_TYPE& out) {
    // Create a float from the binary representation of FPDATA_WORD
    float temp = *reinterpret_cast<const float*>(&in);  // reinterpret cast
    out = FLOAT_TYPE(temp);  // Assign float to FLOAT_TYPE
}

// Function to convert FLOAT_TYPE to FPDATA_WORD
inline void fp2int(const FLOAT_TYPE& in, FPDATA_WORD& out) {
    // Convert FLOAT_TYPE to float and then reinterpret to FPDATA_WORD
    float temp = in.to_float();  // Convert FLOAT_TYPE to float
    out = *reinterpret_cast<FPDATA_WORD*>(&temp);  // reinterpret cast to FPDATA_WORD
}

inline void int2fx(const FPDATA_WORD& in, FPDATA& out)
{ out.set_slc(0,in.slc<FPDATA_WL>(0)); }

inline void fx2int(const FPDATA& in, FPDATA_WORD& out)
{ out.set_slc(0,in.slc<FPDATA_WL>(0)); }

// Convert ac_float to ac_int
inline void float2int(const FLOAT_TYPE& in, FPDATA_WORD& out) {
    // Use the raw bit representation of the FLOAT_TYPE
    out = in.to_ac_int();  // Convert to ac_int directly
}
#endif
