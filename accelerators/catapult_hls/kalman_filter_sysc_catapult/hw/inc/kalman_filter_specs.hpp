// Copyright (c) 2011-2024 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0

#ifndef __ACCSPECS__
#define __ACCSPECS__

#include <mc_connections.h>
#include "nvhls_connections.h"
#include <nvhls_int.h>
#include <nvhls_types.h>
#include <nvhls_vector.h>
#include "kalman_filter_conf_info.hpp"
#include "esp_dma_info_sysc.hpp"
#include <ArbitratedScratchpadDP.h>


// #define iterations 7
// #define dimensions 6

#define STATE_SIZE 6  // Number of states n
#define MEAS_SIZE 164  // Number of measurements m
#define SAMPLES 2  // Number of measurements m

#define iterations SAMPLES
#define dimensions STATE_SIZE
#define const_mat_dim STATE_SIZE


/* <<--defines-->> */
#define DATA_WIDTH 32
#define DMA_SIZE SIZE_WORD
// #define PLM_OUT_WORD 100
// #define PLM_IN_WORD 6400
#define PLM_IN_WORD 32000
#define PLM_OUT_WORD 32000
#define MEM_SIZE 416000/(DMA_WIDTH/8)

#if (DMA_WIDTH == 32)
/* <<--defines_32-->> */
#define DMA_BEAT_PER_WORD 1
#define DMA_WORD_PER_BEAT 1
#define PLM_IN_WP 1
#define PLM_OUT_RP 1
#elif (DMA_WIDTH == 64)
/* <<--defines_64-->> */
#define DMA_BEAT_PER_WORD 1
#define DMA_WORD_PER_BEAT 2
#define PLM_IN_WP 2
#define PLM_OUT_RP 2
#endif

#define PLM_IN_RP 1
#define PLM_OUT_WP 1


const unsigned int inwp = PLM_IN_WP;
const unsigned int inrp = PLM_IN_RP;
const unsigned int outwp = PLM_OUT_WP;
const unsigned int outrp = PLM_OUT_RP;
const unsigned int inbks = PLM_IN_WP;
const unsigned int outbks = PLM_OUT_RP;
const unsigned int inebks = PLM_IN_WORD/PLM_IN_WP;
const unsigned int outebks = PLM_OUT_WORD/PLM_OUT_RP;
const unsigned int in_as = nvhls::index_width<inbks * inebks>::val;
const unsigned int out_as = nvhls::index_width<outbks * outebks>::val;

template<unsigned int kAddressSz>
struct Address{
    typedef NVUINTW(kAddressSz) value;
};

typedef NVUINTW(DATA_WIDTH) DATA_TYPE;

template<unsigned int kAddressSz, unsigned int numports>
class plm_WR : public nvhls_message{
public:
    nvhls::nv_scvector<NVUINTW(kAddressSz), numports> indx;
    nvhls::nv_scvector<DATA_TYPE, numports> data;

    static const unsigned int width = nvhls::nv_scvector<NVUINTW(kAddressSz), numports>::width + nvhls::nv_scvector<DATA_TYPE, numports>::width;
    template <unsigned int Size>
    void Marshall(Marshaller<Size>& m) {
        m & indx;
        m & data;
    }
    plm_WR() {
        Reset();
    }
    void Reset() {
        indx = 0;
        data = 0;
    }
};

template <unsigned int kAddressSz, unsigned int numports>
class plm_RRq : public nvhls_message {
public:
    nvhls::nv_scvector<NVUINTW(kAddressSz), numports> indx;

    static const unsigned int width = nvhls::nv_scvector<NVUINTW(kAddressSz), numports>::width;
    template <unsigned int Size>
    void Marshall(Marshaller<Size> &m) {
        m &indx;
    }

    plm_RRq() {
        Reset();
    }
    void Reset() {
        indx = 0;
    }
};

template <unsigned int numports>
class plm_RRs : public nvhls_message {
public:
    nvhls::nv_scvector<DATA_TYPE, numports> data;
    static const unsigned int width = nvhls::nv_scvector<DATA_TYPE, numports>::width;
    template <unsigned int Size>
    void Marshall(Marshaller<Size> &m) {
        m &data;
    }

    plm_RRs() {
        Reset();
    }
    void Reset() {
        data = 0;
    }
};

#endif
