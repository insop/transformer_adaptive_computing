/*!
 *  Copyright (c) 2018 by Contributors
 * \file vta.h
 * \brief Type definitions and prototype for VTA HLS design.
 */

// ISS based on the file from vta.h and changed to vta_alveo.h

#ifndef VTA_VTA_H_
#define VTA_VTA_H_

#include <ap_axi_sdata.h>
#include <ap_int.h>
#include <assert.h>
#include <hls_stream.h>

#include "./hw_spec.h"

/* \typedef inp_T Input datatype*/
typedef ap_int<VTA_INP_WIDTH> inp_T;

/* \typedef wgt_T Weight datatype*/
typedef ap_int<VTA_WGT_WIDTH> wgt_T;

/* \typedef acc_T Accumulator datatype*/
typedef ap_int<VTA_ACC_WIDTH> acc_T;

/* \typedef mul_T Multiplier output datatype*/
typedef ap_int<VTA_WGT_WIDTH+VTA_INP_WIDTH+1> mul_T;

/* \typedef sum_T GEMM accumulator datatype*/
typedef ap_int<VTA_WGT_WIDTH+VTA_INP_WIDTH+VTA_LOG_BLOCK_IN+1> sum_T;

/* \typedef inp_vec_T Input vector datatype*/
typedef ap_uint<VTA_INP_WIDTH*VTA_BLOCK_IN> inp_vec_T;

/* \typedef wgt_vec_T Weight vector datatype*/
typedef ap_uint<VTA_WGT_WIDTH*VTA_BLOCK_IN> wgt_vec_T;

/* \typedef acc_vec_T Accumulator vector datatype*/
typedef ap_uint<VTA_ACC_WIDTH*VTA_BLOCK_OUT> acc_vec_T;

/* \typedef uop_idx_T Micro-op SRAM index datatype*/
typedef ap_uint<VTA_LOG_UOP_BUFF_DEPTH+1> uop_idx_T;

/* \typedef inp_idx_T Input SRAM index datatype*/
typedef ap_uint<VTA_LOG_INP_BUFF_DEPTH+1> inp_idx_T;

/* \typedef wgt_idx_T Weight SRAM index datatype*/
typedef ap_uint<VTA_LOG_WGT_BUFF_DEPTH+1> wgt_idx_T;

/* \typedef acc_idx_T Accumulator SRAM index datatype*/
typedef ap_uint<VTA_LOG_ACC_BUFF_DEPTH+1> acc_idx_T;

/* \typedef opcode_T Opcode datatype*/
typedef ap_uint<VTA_OPCODE_BIT_WIDTH> opcode_T;

/* \typedef memop_id_T Memory operation ID datatype*/
typedef ap_uint<VTA_MEMOP_ID_BIT_WIDTH> memop_id_T;

/* \typedef memop_sram_T Memory operation SRAM index datatype*/
typedef ap_uint<VTA_MEMOP_SRAM_ADDR_BIT_WIDTH> memop_sram_T;

/* \typedef memop_dram_T Memory operation DRAM index datatype*/
typedef ap_uint<VTA_MEMOP_DRAM_ADDR_BIT_WIDTH> memop_dram_T;

/*!
* \brief Accel design.
* \param insn_count Total instruction count. AXI-lite memory mapped register.
* \param insns Instruction data base address in DRAM. AXI-4 master port.
* \param uops Micro-op data base address in DRAM. AXI-4 master port.
* \param inputs Input data base address in DRAM. AXI-4 master port.
* \param weights Weight data base address in DRAM. AXI-4 master port.
* \param biases Bias data base address in DRAM. AXI-4 master port.
* \param outputs Output data base address in DRAM. AXI-4 master port.
*/
void vta(
  uint32_t insn_count,
  volatile insn_T *insns,
  volatile uop_T *uops,
  volatile inp_vec_T *inputs,
  volatile wgt_vec_T *weights,
  volatile acc_vec_T *biases,
  volatile acc_vec_T *outputs);

#endif  // VTA_VTA_H_
