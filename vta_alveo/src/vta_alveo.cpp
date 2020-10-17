/*!
 *  Copyright (c) 2018 by Contributors
 * \file vta.cc
 * \brief VTA HLS design.
 */

// ISS modefied file from lab1 vta.cc

#include <string.h>

#include "vta_alveo.h"

template <typename WIDE_T, typename NARROW_T, typename IDX_T, int WIDE_W, int NARROW_W, int Y_DIM, int X_DIM>
void read_tensor(
  IDX_T idx,
  WIDE_T src[][NARROW_W * Y_DIM * X_DIM / WIDE_W],
  NARROW_T dst[Y_DIM][X_DIM]) {
#pragma HLS INLINE

  // Read in weight tensor
  for (int p = 0; p < NARROW_W * Y_DIM * X_DIM / WIDE_W; p++) {
    WIDE_T packet = src[idx][p];
    for (int w = 0; w < (WIDE_W / NARROW_W); w++) {
      int x = (p * (WIDE_W / NARROW_W) + w) / X_DIM;
      int y = (p * (WIDE_W / NARROW_W) + w) % X_DIM;
      dst[x][y] = (NARROW_T) packet.range((w + 1) * NARROW_W - 1, w * NARROW_W);
    }
  }
}
template <typename WIDE_T, typename NARROW_T, typename IDX_T, int WIDE_W, int NARROW_W, int Y_DIM, int X_DIM>
void write_tensor(
  IDX_T idx,
  NARROW_T src[Y_DIM][X_DIM],
  WIDE_T dst[][NARROW_W * Y_DIM * X_DIM / WIDE_W]) {
#pragma HLS INLINE

  for (int p = 0; p < NARROW_W * Y_DIM * X_DIM / WIDE_W; p++) {
    WIDE_T packet = 0;
    for (int w = 0; w < (WIDE_W / NARROW_W); w++) {
      int x = (p * (WIDE_W / NARROW_W) + w) / X_DIM;
      int y = (p * (WIDE_W / NARROW_W) + w) % X_DIM;
      packet.range((w + 1) * NARROW_W - 1, w * NARROW_W) = src[x][y];
    }
    dst[idx][p] = packet;
  }
}

extern "C" {

void vta_alveo(
  uint32_t insn_count,
  volatile insn_T *insns,
  volatile uop_T *uops,
  volatile inp_vec_T *inputs,
  volatile wgt_vec_T *weights,
  volatile acc_vec_T *biases,
  volatile acc_vec_T *outputs) {
#pragma HLS INTERFACE s_axilite port = insn_count bundle = CONTROL_BUS
#pragma HLS INTERFACE m_axi port = insns offset = slave bundle = ins_port
#pragma HLS INTERFACE m_axi port = uops offset = slave bundle = uop_port
#pragma HLS INTERFACE m_axi port = inputs offset = slave bundle = narrow_port
#pragma HLS INTERFACE m_axi port = weights offset = slave bundle = narrow_port
#pragma HLS INTERFACE m_axi port = biases offset = slave bundle = wide_port
#pragma HLS INTERFACE m_axi port = outputs offset = slave bundle = wide_port
#pragma HLS INTERFACE s_axilite port = return bundle = CONTROL_BUS

  // Local SRAMs
  uop_T uop_mem[VTA_UOP_BUFF_DEPTH];
  inp_vec_T inp_mem[VTA_INP_BUFF_DEPTH][VTA_BATCH];     // inp_vec_T is VTA_BLOCK_IN wide
  wgt_vec_T wgt_mem[VTA_WGT_BUFF_DEPTH][VTA_BLOCK_OUT]; // wgt_vec_T is VTA_BLOCK_IN wide
  acc_vec_T acc_mem[VTA_ACC_BUFF_DEPTH][VTA_BATCH];     // acc_vec_T is VTA_BLOCK_OUT wide
  // TODO Part 1.b: Insert storage pragmas to improve READ_GEMM_UOP II

  INSN_LOOP: for (int pc = 0; pc < insn_count; pc++) {
    // Load instruction fields
    VTAInsn insn;
    insn.generic = (insn_T) insns[pc];
    // Do some partial decoding
    opcode_T opcode = insn.mem.opcode;
    // Perform appropriate hardware action
    if (opcode == VTA_OPCODE_LOAD || opcode == VTA_OPCODE_STORE) {
      // Decode load instruction
      memop_id_T memory_type = insn.mem.memory_type;
      memop_sram_T sram_base = insn.mem.sram_base;
      memop_dram_T dram_base = insn.mem.dram_base;
      memop_sram_T size = insn.mem.size;
      if (opcode == VTA_OPCODE_LOAD) {
        // Copy the data to the target SRAM
        if (memory_type == VTA_MEM_ID_UOP) {
          memcpy(&uop_mem[sram_base],
                 (const uop_T*) &uops[dram_base],
                 size * sizeof(uop_T));
        } else if (memory_type == VTA_MEM_ID_INP) {
          memcpy(&inp_mem[sram_base][0],
                 (const inp_vec_T*) &inputs[dram_base * VTA_BATCH],
                 size * sizeof(inp_vec_T) * VTA_BATCH);
        } else if (memory_type == VTA_MEM_ID_WGT) {
          memcpy(&wgt_mem[sram_base][0],
                 (const wgt_vec_T*) &weights[dram_base * VTA_BLOCK_OUT],
                 size * sizeof(wgt_vec_T) * VTA_BLOCK_OUT);
        } else if (memory_type == VTA_MEM_ID_ACC) {
          memcpy(&acc_mem[sram_base][0],
                 (const acc_vec_T*) &biases[dram_base * VTA_BATCH],
                 size * sizeof(acc_vec_T) * VTA_BATCH);
        }
      } else {
        // Copy the accumulator data to DRAM
        memcpy(const_cast<acc_vec_T*>(&outputs[dram_base * VTA_BATCH]),
               (const acc_vec_T*) &acc_mem[sram_base][0],
               size * sizeof(acc_vec_T) * VTA_BATCH);
      }
    } else if (opcode == VTA_OPCODE_GEMM) {
      // Decode GEMM instruction
      bool reset_acc = insn.gemm.reset_reg;
      uop_idx_T uop_bgn = insn.gemm.uop_bgn;
      uop_idx_T uop_end = insn.gemm.uop_end;
      // Iterate over micro ops
      READ_GEMM_UOP: for (int upc = uop_bgn; upc < uop_end; upc++) {
        // TODO Part 1.1a: Implement GEMM core
	VTAUop  uop;
	uop.generic = (uop_T) uops[upc];

        acc_idx_T acc_idx = uop.gemm.dst_idx;
        inp_idx_T inp_idx = uop.gemm.src_idx;
        wgt_idx_T wgt_idx = uop.gemm.wgt_idx;

        // Read in weight tensor
        wgt_T w_tensor[VTA_BLOCK_OUT][VTA_BLOCK_IN];
        read_tensor<wgt_vec_T, wgt_T, wgt_idx_T, VTA_WGT_WIDTH*VTA_BLOCK_IN, VTA_WGT_WIDTH, VTA_BLOCK_OUT, VTA_BLOCK_IN>(wgt_idx, wgt_mem, w_tensor);

        // Read in input tensor
        inp_T i_tensor[VTA_BATCH][VTA_BLOCK_IN];
        read_tensor<inp_vec_T, inp_T, inp_idx_T, VTA_INP_WIDTH*VTA_BLOCK_IN, VTA_INP_WIDTH, VTA_BATCH, VTA_BLOCK_IN>(inp_idx, inp_mem, i_tensor);

        // Read in acc tensor
        acc_T a_tensor[VTA_BATCH][VTA_BLOCK_OUT];
        read_tensor<acc_vec_T, acc_T, acc_idx_T, VTA_ACC_WIDTH*VTA_BLOCK_OUT, VTA_ACC_WIDTH, VTA_BATCH, VTA_BLOCK_OUT>(acc_idx, acc_mem, a_tensor);

	for (int i = 0; i < VTA_BATCH; i++) {
	  for(int j = 0; j < VTA_BLOCK_OUT; j++) {

            acc_T accum = a_tensor[i][j];
            acc_T tmp = 0;
	    for (int k = 0; k < VTA_BLOCK_IN; k++) {
              wgt_T w_elem = w_tensor[j][k];
              inp_T i_elem = i_tensor[i][k];
              acc_T prod = i_elem * w_elem;
              tmp += (acc_T) prod;
            }
	    accum += tmp;
	    a_tensor[i][j] = reset_acc == true ? (acc_T)0 : accum;
	  }
	}
        // Write the results back into accumulator
        write_tensor<acc_vec_T, acc_T, acc_idx_T, VTA_ACC_WIDTH*VTA_BLOCK_OUT, VTA_ACC_WIDTH, VTA_BATCH, VTA_BLOCK_OUT>(acc_idx, a_tensor, acc_mem);
	
        // TODO Part 1.1b: Insert loop pragmas to improve READ_GEMM_UOP II
#pragma HLS PIPELINE II=1
#pragma HLS dependence variable=w_tensor inter false
#pragma HLS dependence variable=i_tensor inter false
#pragma HLS dependence variable=a_tensor inter false
      }
    }
  }
}

} // extern "C"
