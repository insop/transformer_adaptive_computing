/*!
 *  Copyright (c) 2018 by Contributors
 * \file hw_spec.h
 * \brief Preprocessor definitions for VTA HLS design (Xilinx independent).
 */

#ifndef VTA_HW_SPEC_H_
#define VTA_HW_SPEC_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

/*! log2 of input type width */
#define VTA_LOG_INP_WIDTH 3
/*! log2 of weight type width */
#define VTA_LOG_WGT_WIDTH 3
/*! log2 of accum type width */
#define VTA_LOG_ACC_WIDTH 5
/*! log2 of batch tensorization size (A in (A,B)x(B,C) mat mul) */
#define VTA_LOG_BATCH 0
/*! log2 of inner block tensorization size (B in (A,B)x(B,C) mat mul) */
#define VTA_LOG_BLOCK_IN 4
/*! log2 of outer block tensorization size (C in (A,B)x(B,C) mat mul) */
#define VTA_LOG_BLOCK_OUT 4
/*! log2 of on-chip micro-op buffer size in B */
#define VTA_LOG_UOP_BUFF_SIZE 14
/*! log2 of on-chip input buffer size in B */
#define VTA_LOG_INP_BUFF_SIZE 15
/*! log2 of on-chip weigth buffer size in B */
#define VTA_LOG_WGT_BUFF_SIZE 18
/*! log2 of on-chip accum buffer size in B */
#define VTA_LOG_ACC_BUFF_SIZE 17

/*! log2 of instruction data type width */
#define VTA_LOG_INS_WIDTH 6
/*! Instruction data type width */
#define VTA_INS_WIDTH (1 << VTA_LOG_INS_WIDTH)
/*! log2 of micro op data type width */
#define VTA_LOG_UOP_WIDTH 5
/*! Micro Op data type width */
#define VTA_UOP_WIDTH (1 << VTA_LOG_UOP_WIDTH)
/*! Weight data type width */
#define VTA_WGT_WIDTH (1 << VTA_LOG_WGT_WIDTH)
/*! Input data type width */
#define VTA_INP_WIDTH (1 << VTA_LOG_INP_WIDTH)
/*! Accumulator data type width */
#define VTA_ACC_WIDTH (1 << VTA_LOG_ACC_WIDTH)

/*! Batch size (corresponds to A in (A,B)x(B,C) mat mult)*/
#define VTA_BATCH (1 << VTA_LOG_BATCH)
/*! Blocking factor of inner most loop (corresponds to B in (A,B)x(B,C) mat mult) */
#define VTA_BLOCK_IN (1 << VTA_LOG_BLOCK_IN)
/*! Blocking factor of the outer loop (corresponds to C in (A,B)x(B,C) mat mult) */
#define VTA_BLOCK_OUT (1 << VTA_LOG_BLOCK_OUT)

/*! Weight vector width */
#define VTA_WGT_VECTOR_WIDTH (VTA_WGT_WIDTH * VTA_BLOCK_IN)
/*! Input vector width */
#define VTA_INP_VECTOR_WIDTH (VTA_INP_WIDTH * VTA_BLOCK_IN)
/*! Accumulator vector width */
#define VTA_ACC_VECTOR_WIDTH (VTA_ACC_WIDTH * VTA_BLOCK_OUT)

/*! On-chip micro-op buffer size in B */
#define VTA_UOP_BUFF_SIZE (1 << VTA_LOG_UOP_BUFF_SIZE)
/*! On-chip weight buffer size in B */
#define VTA_WGT_BUFF_SIZE (1 << VTA_LOG_WGT_BUFF_SIZE)
/*! On-chip activation buffer size in B */
#define VTA_INP_BUFF_SIZE (1 << VTA_LOG_INP_BUFF_SIZE)
/*! On-chip accumulator buffer size in B */
#define VTA_ACC_BUFF_SIZE (1 << VTA_LOG_ACC_BUFF_SIZE)

/*! Size of instruction buffer element in B */
#define VTA_INS_ELEM_BYTES (VTA_INS_WIDTH / 8)
/*! Size of uop buffer element in B*/
#define VTA_UOP_ELEM_BYTES (VTA_UOP_WIDTH / 8)
/*! Size of activation buffer element in B*/
#define VTA_INP_ELEM_BYTES (VTA_BATCH * VTA_BLOCK_IN * VTA_INP_WIDTH / 8)
/*! Size of weight buffer element in B*/
#define VTA_WGT_ELEM_BYTES (VTA_BLOCK_OUT * VTA_BLOCK_IN * VTA_WGT_WIDTH / 8)
/*! Size of accumulator buffer element in B*/
#define VTA_ACC_ELEM_BYTES (VTA_BATCH * VTA_BLOCK_OUT * VTA_ACC_WIDTH / 8)

/*! On-chip micro-op buffer depth */
#define VTA_UOP_BUFF_DEPTH (VTA_UOP_BUFF_SIZE / VTA_UOP_ELEM_BYTES)
/*! log2 of on-chip micro-op buffer depth */
#define VTA_LOG_UOP_BUFF_DEPTH (VTA_LOG_UOP_BUFF_SIZE - VTA_LOG_UOP_WIDTH + 3)
// ! \brief On-chip weight buffer depth
#define VTA_WGT_BUFF_DEPTH (VTA_WGT_BUFF_SIZE / VTA_WGT_ELEM_BYTES)
/*! log2 of weight micro-op buffer depth */
#define VTA_LOG_WGT_BUFF_DEPTH \
    (VTA_LOG_WGT_BUFF_SIZE - VTA_LOG_BLOCK_OUT - VTA_LOG_BLOCK_IN - VTA_LOG_WGT_WIDTH + 3)
/*! On-chip activation buffer depth */
#define VTA_INP_BUFF_DEPTH (VTA_INP_BUFF_SIZE / VTA_INP_ELEM_BYTES)
/*! log2 of activation micro-op buffer depth */
#define VTA_LOG_INP_BUFF_DEPTH \
    (VTA_LOG_INP_BUFF_SIZE - VTA_LOG_BATCH - VTA_LOG_BLOCK_IN - VTA_LOG_INP_WIDTH + 3)
/*! On-chip accumulator buffer depth */
#define VTA_ACC_BUFF_DEPTH (VTA_ACC_BUFF_SIZE / VTA_ACC_ELEM_BYTES)
/*! log2 of on-chip accumulator buffer depth */
#define VTA_LOG_ACC_BUFF_DEPTH \
    (VTA_LOG_ACC_BUFF_SIZE - VTA_LOG_BATCH - VTA_LOG_BLOCK_OUT - VTA_LOG_ACC_WIDTH + 3)

/*! Instruction opcode field bitwidth */
#define VTA_OPCODE_BIT_WIDTH 2

/*! Opcode: load encoding */
#define VTA_OPCODE_LOAD 0
/*! Opcode: store encoding */
#define VTA_OPCODE_STORE 1
/*! Opcode: GEMM encoding */
#define VTA_OPCODE_GEMM 2

/*! Memory type field bitwidth */
#define VTA_MEMOP_ID_BIT_WIDTH 2
/*! Load/Store Instruction: DRAM address width*/
#define VTA_MEMOP_SRAM_ADDR_BIT_WIDTH 14
/*! Load/Store Instruction: DRAM address width*/
#define VTA_MEMOP_DRAM_ADDR_BIT_WIDTH 32

/*! Mem ID constant: uop memory */
#define VTA_MEM_ID_UOP 0
/*! Mem ID constant: weight memory */
#define VTA_MEM_ID_WGT 1
/*! Mem ID constant: input memory */
#define VTA_MEM_ID_INP 2
/*! Mem ID constant: accumulator/bias memory */
#define VTA_MEM_ID_ACC 3

// Instruction organization layout:
//
// LOAD/STORE
// _____________________________|_type______________|
// arg 0: opcode                | opcode_T          |
// arg 1: memory_type           | memop_id_T        |
// arg 2: sram_base             | memop_sram_T      |
// arg 3: dram_base             | memop_dram_T      |
// arg 4: size                  | memop_size_T      |
//
// GEMM
// _____________________________|_type______________|
// arg 0: opcode                | opcode_T          |
// arg 1: reset_reg             | bool              |
// arg 2: uop_bgn               | uop_idx_T         |
// arg 3: uop_end               | uop_idx_T         |

typedef uint64_t VTAGenericInsn;

/*! \brief VTA load/store instruction */
typedef struct {
  /*! \brief The instruction opcode */
  uint64_t opcode         : VTA_OPCODE_BIT_WIDTH;
  /*! \brief Source/destination SRAM for store/load instruction */
  uint64_t memory_type    : VTA_MEMOP_ID_BIT_WIDTH;
  /*! \brief SRAM base address (pointer to memory elem type) */
  uint64_t sram_base      : VTA_MEMOP_SRAM_ADDR_BIT_WIDTH;
  /*! \brief DRAM base address (pointer to memory elem type) */
  uint64_t dram_base      : VTA_MEMOP_DRAM_ADDR_BIT_WIDTH;
  /*! \brief Transfer size (in elems) */
  uint64_t size           : VTA_MEMOP_SRAM_ADDR_BIT_WIDTH;
} VTAMemInsn;

/*! \brief VTA GEMM instruction */
typedef struct {
  /*! \brief The instruction opcode */
  uint64_t opcode         : VTA_OPCODE_BIT_WIDTH;
  /*! \brief Reset register */
  uint64_t reset_reg      : 1;
  /*! \brief Micro-op begin address */
  uint64_t uop_bgn        : VTA_LOG_UOP_BUFF_DEPTH;
  /*! \brief Micro-op end address */
  uint64_t uop_end        : VTA_LOG_UOP_BUFF_DEPTH + 1;
} VTAGemInsn;

/*! \brief VTA ALU instruction converter */
union VTAInsn {
  /*! \brief VTA generic instruction */
  VTAGenericInsn generic;
  /*! \brief VTA load/store instruction */
  VTAMemInsn mem;
  /*! \brief VTA GEMM instruction */
  VTAGemInsn gemm;
};

typedef uint32_t VTAGenericUop;

/*! \brief VTA micro-op for GEMM/ALU instruction */
typedef struct {
  /*! \brief Destination index (indexes accum buffer) */
  uint32_t dst_idx    : VTA_LOG_ACC_BUFF_DEPTH;
  /*! \brief Source index (indexes input buffer for GEMM or accum buffer for ALU) */
  uint32_t src_idx    : VTA_LOG_INP_BUFF_DEPTH;
  /*! \brief Weight index (indexes weight buffer) */
  uint32_t wgt_idx    : VTA_LOG_WGT_BUFF_DEPTH;
} VTAGemUop;


/*! \brief VTA ALU micro-op converter */
union VTAUop {
  /*! \brief VTA generic micro-op */
  VTAGenericUop generic;
  /*! \brief VTA GEMM micro-op */
  VTAGemUop gemm;
};

/* \typedef insn_T Instruction datatype*/
typedef VTAGenericInsn insn_T;
/* \typedef uop_T Micro-op datatype*/
typedef VTAGenericUop uop_T;

#ifdef __cplusplus
}
#endif
#endif  // VTA_HW_SPEC_H_
