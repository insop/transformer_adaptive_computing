/*!
 *  Copyright (c) 2018 by Contributors
 * \file test_lib.cpp
 * \brief Test library for the VTA design simulation and driver tests.
 */

#include "./test_lib.h"

uint32_t globalSeed;

#define VTA_DEBUG 1

#ifdef NO_SIM

// pync driver to call h/w
uint64_t vta_alveo(uint32_t insn_count,
             VTAGenericInsn *insns,
             VTAGenericUop *uops,
             inp_T *inputs,
             wgt_T *weights,
             acc_T *biases,
             acc_T *outputs) {
  // Performance counter variables
  uint64_t t_fpga;
  struct timespec start, stop;

#if 0
  // Program
  Program("vta.bit");

  // HLS Design Handle
  void* hls_handle = MapRegister(HLS_IP_ADDR, RANGE);

  // Physical address pointers
  uint32_t insn_phy = insns ? cma_get_phy_addr(insns) : 0;
  uint32_t uop_phy = uops ? cma_get_phy_addr(uops) : 0;
  uint32_t input_phy = inputs ? cma_get_phy_addr(inputs) : 0;
  uint32_t weight_phy = weights ? cma_get_phy_addr(weights) : 0;
  uint32_t bias_phy = biases ? cma_get_phy_addr(biases) : 0;
  uint32_t output_phy = outputs ? cma_get_phy_addr(outputs) : 0;
  clock_gettime(CLOCK_REALTIME, &start);

  WriteMappedReg(hls_handle, 0x10, insn_count);
  WriteMappedReg(hls_handle, 0x18, insn_phy);
  WriteMappedReg(hls_handle, 0x20, uop_phy);
  WriteMappedReg(hls_handle, 0x28, input_phy);
  WriteMappedReg(hls_handle, 0x30, weight_phy);
  WriteMappedReg(hls_handle, 0x38, bias_phy);
  WriteMappedReg(hls_handle, 0x40, output_phy);
  WriteMappedReg(hls_handle, 0x00, START);

  // Loop until done
  unsigned t, flag = 0;
  unsigned wait_cycles = 1000000;
  for (t = 0; t < wait_cycles; ++t) {
    flag = ReadMappedReg(hls_handle, 0x00);
    if (flag & DONE) break;
  }

  clock_gettime(CLOCK_REALTIME, &stop);
  t_fpga = 1000000000ULL * (stop.tv_sec - start.tv_sec) + (stop.tv_nsec - start.tv_nsec);
#endif

  return t_fpga;
}

#endif  // NO_SIM

template <typename T, int T_WIDTH>
void pack2dBuffer(T *dst, T **src, int y_size, int x_size, int y_block, int x_block) {
  int buffer_idx = 0;
  for (int i = 0; i < y_size / y_block; i++) {
    for (int j = 0; j < x_size / x_block; j++) {
      for (int k = 0; k < y_block; k++) {
        if (T_WIDTH < 8) {
          for (int l = 0; l < x_block; l += 8 / T_WIDTH) {
            dst[buffer_idx] = 0;
            for (int m = 0; m < 8 / T_WIDTH; m++) {
              dst[buffer_idx] |= (src[i * y_block + k][j * x_block + l + m] &
                ((1ULL << T_WIDTH) - 1)) << (m * T_WIDTH);
            }
            buffer_idx++;
          }
        } else {
          for (int l = 0; l < x_block; l++) {
            dst[buffer_idx++] = src[i * y_block + k][j * x_block + l];
          }
        }
      }
    }
  }
}

template <typename T, int T_WIDTH>
void unpack2dBuffer(T **dst, T *src, int y_size, int x_size, int y_block, int x_block) {
  int buffer_idx = 0;
  for (int i = 0; i < y_size / y_block; i++) {
    for (int j = 0; j < x_size / x_block; j++) {
      for (int k = 0; k < y_block; k++) {
        if (T_WIDTH < 8) {
          for (int l = 0; l < x_block; l += 8 / T_WIDTH) {
            for (int m = 0; m < 8 / T_WIDTH; m++) {
              dst[i * y_block + k][j * x_block + l + m] = (int)
                (((long long int) src[buffer_idx]) >> (m * T_WIDTH))
                & ((1ULL << T_WIDTH) - 1);
            }
            buffer_idx++;
          }
        } else {
          for (int l = 0; l < x_block; l++) {
            dst[i * y_block + k][j * x_block + l] = src[buffer_idx++];
          }
        }
      }
    }
  }
}

template <typename T, int T_WIDTH>
void pack4dBuffer(T *dst, T ****src, int ax0, int ax1, int ax2, int ax3,
    int ax0_block, int ax3_block) {
  int buffer_idx = 0;
  for (int i = 0; i < ax0 / ax0_block; i++) {
    for (int j = 0; j < ax1; j++) {
      for (int k = 0; k < ax2; k++) {
        for (int l = 0; l < ax3 / ax3_block; l++) {
          for (int di = 0; di < ax0_block; di++) {
            if (T_WIDTH < 8) {
              for (int dl = 0; dl < ax3_block; dl += 8 / T_WIDTH) {
                dst[buffer_idx] = 0;
                for (int m = 0; m < 8 / T_WIDTH; m++) {
                  dst[buffer_idx] |= (src[i * ax0_block + di][j][k][l * ax3_block + dl + m] &
                    ((1ULL << T_WIDTH) - 1)) << (m * T_WIDTH);
                }
                buffer_idx++;
              }
            } else {
              for (int dl = 0; dl < ax3_block; dl++) {
                dst[buffer_idx++] = src[i * ax0_block + di][j][k][l * ax3_block + dl];
              }
            }
          }
        }
      }
    }
  }
}

template <typename T, int T_WIDTH>
void unpack4dBuffer(T ****dst, T *src, int ax0, int ax1, int ax2, int ax3,
    int ax0_block, int ax3_block) {
  int buffer_idx = 0;
  for (int i = 0; i < ax0 / ax0_block; i++) {
    for (int j = 0; j < ax1; j++) {
      for (int k = 0; k < ax2; k++) {
        for (int l = 0; l < ax3 / ax3_block; l++) {
          for (int di = 0; di < ax0_block; di++) {
            if (T_WIDTH < 8) {
              for (int dl = 0; dl < ax3_block; dl += 8 / T_WIDTH) {
                for (int m = 0; m < 8 / T_WIDTH; m++) {
                  dst[i * ax0_block + di][j][k][l * ax3_block + dl + m] = (int)
                    (((long long int) src[buffer_idx]) >> (m * T_WIDTH))
                    & ((1ULL << T_WIDTH) - 1);
                }
                buffer_idx++;
              }
            } else {
              for (int dl = 0; dl < ax3_block; dl++) {
                dst[i * ax0_block + di][j][k][l * ax3_block + dl] = src[buffer_idx++];
              }
            }
          }
        }
      }
    }
  }
}

template <typename T, int T_WIDTH>
T ** allocInit2dArray(int rows, int cols) {
  // Allocate
  T **array = static_cast<T **>(malloc(sizeof(T *) * rows));
  for (int i = 0; i < rows; i++) {
    array[i] = static_cast<T *>(malloc(sizeof(T) * cols));
  }
  // Init
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      array[i][j] =
          static_cast<T>(rand_r(&globalSeed) % (1LL << (T_WIDTH - 1)) - (1LL << (T_WIDTH - 2)));
    }
  }
  return array;
}

template <typename T>
T ** alloc2dArray(int rows, int cols) {
  T **array = static_cast<T **>(malloc(sizeof(T *) * rows));
  for (int i = 0; i < rows; i++) {
    array[i] = static_cast<T *>(malloc(sizeof(T) * cols));
  }
  return array;
}

template <typename T>
void free2dArray(T **array, int rows, int cols) {
  for (int i = 0; i < rows; i++) {
    free(array[i]);
  }
  free(array);
}

template <typename T, int T_WIDTH>
T **** allocInit4dArray(int ax0, int ax1, int ax2, int ax3) {
  // Allocate
  T ****array = static_cast<T ****>(malloc(sizeof(T ***) * ax0));
  for (int i = 0; i < ax0; i++) {
    array[i] = static_cast<T ***>(malloc(sizeof(T **) * ax1));
    for (int j = 0; j < ax1; j++) {
      array[i][j] = static_cast<T **>(malloc(sizeof(T *) * ax2));
      for (int k = 0; k < ax2; k++) {
        array[i][j][k] = static_cast<T *>(malloc(sizeof(T) * ax3));
      }
    }
  }
  // Init
  for (int i = 0; i < ax0; i++) {
    for (int j = 0; j < ax1; j++) {
      for (int k = 0; k < ax2; k++) {
        for (int l = 0; l < ax3; l++) {
          array[i][j][k][l] =
              static_cast<T>(rand_r(&globalSeed) % (1LL << (T_WIDTH - 1)) - (1LL << (T_WIDTH - 2)));
        }
      }
    }
  }
  return array;
}

template <typename T>
T **** alloc4dArray(int ax0, int ax1, int ax2, int ax3) {
  T ****array = static_cast<T ****>(malloc(sizeof(T ***) * ax0));
  for (int i = 0; i < ax0; i++) {
    array[i] = static_cast<T ***>(malloc(sizeof(T **) * ax1));
    for (int j = 0; j < ax1; j++) {
      array[i][j] = static_cast<T **>(malloc(sizeof(T *) * ax2));
      for (int k = 0; k < ax2; k++) {
        array[i][j][k] = static_cast<T *>(malloc(sizeof(T) * ax3));
      }
    }
  }
  return array;
}

template <typename T>
void free4dArray(T ****array, int ax0, int ax1, int ax2, int ax3) {
  for (int i = 0; i < ax0; i++) {
    for (int j = 0; j < ax1; j++) {
      for (int k = 0; k < ax2; k++) {
        free(array[i][j][k]);
      }
      free(array[i][j]);
    }
    free(array[i]);
  }
  free(array);
}

void * allocBuffer(size_t num_bytes) {
#ifdef NO_SIM
  // ISS
  return malloc(num_bytes);
  //return cma_alloc(num_bytes, CACHED);
#else
  return malloc(num_bytes);
#endif
}

void freeBuffer(void * buffer) {
#ifdef NO_SIM
  // ISS
  return free(buffer);
  //return cma_free(buffer);
#else
  return free(buffer);
#endif
}

VTAGenericInsn getLoadStoreInsn(int opcode, int type, int size) {
  // Converter
  union VTAInsn converter;
  // Memory instruction initialization
  VTAMemInsn insn = {};
  insn.opcode = opcode;
  insn.memory_type = type;
  insn.sram_base = 0;
  insn.dram_base = 0;
  insn.size = size;
  converter.mem = insn;
  return converter.generic;
}

VTAGenericInsn getGEMMInsn(int uop_size, bool reset) {
  // Converter
  union VTAInsn converter;
  // GEVM instruction initialization
  VTAGemInsn insn;
  insn.opcode = VTA_OPCODE_GEMM;
  insn.reset_reg = reset;
  insn.uop_bgn = 0;
  insn.uop_end = uop_size;
  converter.gemm = insn;
  return converter.generic;
}

VTAGenericUop * getResetUops(int batch, int out_feat) {
  // Converter
  union VTAUop converter;
  // Derive the total uop size
  int uop_size = batch * out_feat;
  // Allocate buffer
#ifdef NO_SIM
  // ISS
  VTAGenericUop *uop_buf = static_cast<VTAGenericUop *>(malloc(sizeof(VTAGenericUop) * uop_size));
  //VTAGenericUop *uop_buf = static_cast<VTAGenericUop *>(cma_alloc(sizeof(VTAGenericUop) * uop_size, CACHED));
#else
  VTAGenericUop *uop_buf = static_cast<VTAGenericUop *>(malloc(sizeof(VTAGenericUop) * uop_size));
#endif
  // Generate micro ops
  int uop_idx = 0;
  for (int i = 0; i < batch; i++) {
    for (int j = 0; j < out_feat; j++) {
      converter.gemm.dst_idx = i * out_feat + j;
      converter.gemm.src_idx = 0;
      converter.gemm.wgt_idx = 0;
      uop_buf[uop_idx++] = converter.generic;
    }
  }
  return uop_buf;
}

VTAGenericUop * getGEMMUops(int batch, int in_feat, int out_feat) {
  // Converter
  union VTAUop converter;
  // Derive the total uop size
  int uop_size = batch * in_feat * out_feat;
  // Allocate buffer
#ifdef NO_SIM
  // ISS
  VTAGenericUop *uop_buf = static_cast<VTAGenericUop *>(malloc(sizeof(VTAGenericUop) * uop_size));
  //VTAGenericUop *uop_buf = static_cast<VTAGenericUop *>(cma_alloc(sizeof(VTAGenericUop) * uop_size, CACHED));
#else
  VTAGenericUop *uop_buf = static_cast<VTAGenericUop *>(malloc(sizeof(VTAGenericUop) * uop_size));
#endif
  // Generate micro ops
  int uop_idx = 0;
  for (int i = 0; i < batch; i++) {
    for (int j = 0; j < in_feat; j++) {
      for (int k = 0; k < out_feat; k++) {
        converter.gemm.dst_idx = i * out_feat + k;
        converter.gemm.src_idx = i * in_feat + j;
        converter.gemm.wgt_idx = k * in_feat + j;
        uop_buf[uop_idx++] = converter.generic;
      }
    }
  }
  return uop_buf;
}

VTAGenericUop * getConv2dUops(int batch, int height, int width, int kheight, int kwidth,
    int in_channels, int out_channels) {
  // Converter
  union VTAUop converter;
  // Derive the total uop size
  int uop_size = batch * in_channels * out_channels * height * width * kheight * kwidth;
  // Derive padding
  int ypad = (kheight - 1) / 2;
  int xpad = (kwidth - 1) / 2;
  // Allocate buffer
#ifdef NO_SIM
   // ISS
  VTAGenericUop *uop_buf = static_cast<VTAGenericUop *>(malloc(sizeof(VTAGenericUop) * uop_size));
  // VTAGenericUop *uop_buf = static_cast<VTAGenericUop *>(cma_alloc(sizeof(VTAGenericUop) * uop_size, CACHED));
#else
  VTAGenericUop *uop_buf = static_cast<VTAGenericUop *>(malloc(sizeof(VTAGenericUop) * uop_size));
#endif
  // Generate micro ops
  int uop_idx = 0;
  for (int b = 0; b < batch; b++) {
    for (int y = ypad; y < height - ypad; y++) {
      for (int x = xpad; x < width - xpad; x++) {
        for (int oc = 0; oc < out_channels; oc++) {
          for (int dy = 0 - ypad; dy < ypad + 1; dy++) {
            for (int dx = 0 - xpad; dx < xpad + 1; dx++) {
              for (int ic = 0; ic < in_channels; ic++) {
                // TODO Part 2: Derive 2D convolution indices
                converter.gemm.dst_idx = 0;
                converter.gemm.src_idx = 0;
                converter.gemm.wgt_idx = 0;
                uop_buf[uop_idx++] = converter.generic;
              }
            }
          }
        }
      }
    }
  }
  return uop_buf;
}

void printParameters() {
  // Some debugging code
  printf("Size of VTAInsn: %d\n", sizeof(VTAGenericInsn));
  printf("Size of VTAUop: %d\n", sizeof(VTAGenericUop));
  printf("VTA_UOP_BUFF_DEPTH: %d\n", VTA_UOP_BUFF_DEPTH);
  printf("VTA_LOG_UOP_BUFF_DEPTH: %d\n", VTA_LOG_UOP_BUFF_DEPTH);
  printf("VTA_WGT_BUFF_DEPTH: %d\n", VTA_WGT_BUFF_DEPTH);
  printf("VTA_LOG_WGT_BUFF_DEPTH: %d\n", VTA_LOG_WGT_BUFF_DEPTH);
  printf("VTA_INP_BUFF_DEPTH: %d\n", VTA_INP_BUFF_DEPTH);
  printf("VTA_LOG_INP_BUFF_DEPTH: %d\n", VTA_LOG_INP_BUFF_DEPTH);
  printf("VTA_ACC_BUFF_DEPTH: %d\n", VTA_ACC_BUFF_DEPTH);
  printf("VTA_LOG_ACC_BUFF_DEPTH: %d\n", VTA_LOG_ACC_BUFF_DEPTH);
  printf("VTA_WGT_WORDS: %d\n", VTA_WGT_BUFF_DEPTH*VTA_BLOCK_IN*VTA_BLOCK_OUT);
  printf("VTA_INP_WORDS: %d\n", VTA_INP_BUFF_DEPTH*VTA_BLOCK_IN);
  printf("VTA_ACC_WORDS: %d\n", VTA_ACC_BUFF_DEPTH*VTA_BLOCK_OUT);
  printf("VTA_INS_ELEM_BYTES: %d\n", VTA_INS_ELEM_BYTES);
  printf("VTA_UOP_ELEM_BYTES: %d\n", VTA_UOP_ELEM_BYTES);
  printf("VTA_INP_ELEM_BYTES: %d\n", VTA_INP_ELEM_BYTES);
  printf("VTA_WGT_ELEM_BYTES: %d\n", VTA_WGT_ELEM_BYTES);
  printf("VTA_ACC_ELEM_BYTES: %d\n", VTA_ACC_ELEM_BYTES);
  printf("VTA_BLOCK_IN: %d\n", VTA_BLOCK_IN);
  printf("VTA_BLOCK_OUT: %d\n", VTA_BLOCK_OUT);
}

void printInstruction(int num_insn, VTAGenericInsn *insns) {
  // Converter
  union VTAInsn c;
  // Iterate over all instructions
  printf("DEBUG - There are %u instructions\n", num_insn);
  for (int i = 0; i < num_insn; i++) {
    // Fetch instruction and decode opcode
    c.generic = insns[i];
    printf("DEBUG - INSTRUCTION %u: ", i);
    if (c.mem.opcode == VTA_OPCODE_LOAD || c.mem.opcode == VTA_OPCODE_STORE) {
      // Print instruction field information
      if (c.mem.opcode == VTA_OPCODE_LOAD) {
        printf("LOAD ");
        if (c.mem.memory_type == VTA_MEM_ID_UOP) printf("UOP\n");
        if (c.mem.memory_type == VTA_MEM_ID_WGT) printf("WGT\n");
        if (c.mem.memory_type == VTA_MEM_ID_INP) printf("INP\n");
        if (c.mem.memory_type == VTA_MEM_ID_ACC) printf("ACC\n");
      }
      if (c.mem.opcode == VTA_OPCODE_STORE) {
        printf("STORE ACC\n");
      }
      printf("\tDRAM: 0x%08x, SRAM:0x%04x\n",
             static_cast<int>(c.mem.dram_base),
             static_cast<int>(c.mem.sram_base));
      printf("\tsize=%d\n", static_cast<int>(c.mem.size));
    } else if (c.mem.opcode == VTA_OPCODE_GEMM) {
      // Print instruction field information
      printf("GEVM\n");
      printf("\treset_out: %d\n", static_cast<int>(c.gemm.reset_reg));
      printf("\trange (%d, %d)\n",
             static_cast<int>(c.gemm.uop_bgn),
             static_cast<int>(c.gemm.uop_end));
    }
  }
}

// Helper function: Print micro-ops status
void printMicroOp(int num_uop, VTAGenericUop *uops) {
  // Converter
  union VTAUop c;
  // Iterate over all micro ops
  printf("DEBUG - There are %u micro-ops\n", num_uop);
  for (int i = 0; i < num_uop; i++) {
    // Read micro-op
    c.generic = uops[i];
    printf("DEBUG - UOP %u: ", i);
    printf("acc=%u, inp=%u, wgt=%u\n", c.gemm.dst_idx, c.gemm.src_idx, c.gemm.wgt_idx);
  }
}

int mem_test(int batch, int out_channels) {
  // Some assertions
  assert(out_channels % VTA_BLOCK_OUT == 0);
  assert(batch % VTA_BATCH == 0);

  printf("=====================================================================================\n");
  printf("INFO - Load/Store test: batch=%d, out_channels=%d\n",
         batch, out_channels);

  // Derive number of elements that need to be loaded/stored
  int ins_size = 2;
  int xfer_size = batch / VTA_BATCH * out_channels / VTA_BLOCK_OUT;
  // Make sure we don't exceed buffer bounds
  assert(xfer_size <= VTA_ACC_BUFF_DEPTH);

  // Initialize instruction buffer
  VTAGenericInsn *insn_buf =
      static_cast<VTAGenericInsn *>(allocBuffer(sizeof(VTAGenericInsn) * ins_size));

  // Load data block
  insn_buf[0] = getLoadStoreInsn(VTA_OPCODE_LOAD,
                                 VTA_MEM_ID_ACC,
                                 xfer_size);
  // Store data block
  insn_buf[1] = getLoadStoreInsn(VTA_OPCODE_STORE,
                                 VTA_MEM_ID_ACC,
                                 xfer_size);

#if VTA_DEBUG == 1
  printInstruction(ins_size, insn_buf);
#endif

  // Initialize input data
  acc_T **inputs = allocInit2dArray<acc_T, VTA_ACC_WIDTH>(batch, out_channels);

  // Reference output
  acc_T **outputs_ref = alloc2dArray<acc_T>(batch, out_channels);
  for (int i = 0; i < batch; i++) {
    for (int j = 0; j < out_channels; j++) {
      outputs_ref[i][j] = inputs[i][j];
    }
  }

  // Prepare the input buffer
  acc_T *input_buf = static_cast<acc_T *>(allocBuffer(VTA_ACC_ELEM_BYTES * xfer_size));
  pack2dBuffer<acc_T, VTA_ACC_WIDTH>(input_buf,
                                     inputs,
                                     batch,
                                     out_channels,
                                     VTA_BATCH,
                                     VTA_BLOCK_OUT);
  // Prepare the output buffer
  acc_T *output_buf = static_cast<acc_T *>(allocBuffer(VTA_ACC_ELEM_BYTES * xfer_size));

#ifdef NO_SIM
  // Invoke the VTA
  uint64_t t_fpga = vta_alveo(ins_size,
                        insn_buf,
                        NULL,
                        NULL,
                        NULL,
                        input_buf,
                        output_buf);
  // Report on timining
  printf("INFO - Synchronization time: %.3lfms\n", static_cast<float>(t_fpga) / 1E6);
  printf("INFO - Throughput: %.3lfGbs/s\n",
         static_cast<float>(xfer_size) * 2 * VTA_ACC_ELEM_BYTES * 8 / t_fpga);
#else
  // Invoke the VTA
  vta_alveo(ins_size,
      (volatile insn_T *) insn_buf,
      NULL,
      NULL,
      NULL,
      (volatile acc_vec_T *) input_buf,
      (volatile acc_vec_T *) output_buf);
#endif

  // Unpack output data
  acc_T **outputs = alloc2dArray<acc_T>(batch, out_channels);
  unpack2dBuffer<acc_T, VTA_ACC_WIDTH>(outputs,
                                       output_buf,
                                       batch,
                                       out_channels,
                                       VTA_BATCH,
                                       VTA_BLOCK_OUT);

  // Correctness checks
  int err = 0;
  for (int i = 0; i < batch; i++) {
    for (int j = 0; j < out_channels; j++) {
      if (outputs_ref[i][j] != outputs[i][j]) {
        err++;
#if VTA_DEBUG == 1
        printf("DEBUG - %d, %d: expected 0x%x but got 0x%x\n", i, j,
               static_cast<int>(outputs_ref[i][j]),
               static_cast<int>(outputs[i][j]));
#endif  // VTA_DEBUG
      }
    }
  }

  // Free all allocated arrays
  free2dArray<acc_T>(inputs, batch, out_channels);
  free2dArray<acc_T>(outputs_ref, batch, out_channels);
  free2dArray<acc_T>(outputs, batch, out_channels);
  freeBuffer(insn_buf);
  freeBuffer(input_buf);
  freeBuffer(output_buf);

  if (err == 0) {
    printf("INFO - Load/Store test successful!\n");
    return 0;
  } else {
    printf("INFO - Load/Store test failed, got %d errors!\n", err);
    return -1;
  }
}

int reset_test(int batch, int out_channels) {
  // Some assertions
  assert(out_channels % VTA_BLOCK_OUT == 0);
  assert(batch % VTA_BATCH == 0);

  printf("=====================================================================================\n");
  printf("INFO - Reset test: batch=%d, out_channels=%d\n",
         batch, out_channels);

  // Derive number of elements that need to be loaded/stored
  int ins_size = 4;
  int uop_size = batch / VTA_BATCH * out_channels / VTA_BLOCK_OUT;
  int xfer_size = batch / VTA_BATCH * out_channels / VTA_BLOCK_OUT;
  // Make sure we don't exceed buffer bounds
  assert(uop_size <= VTA_UOP_BUFF_DEPTH);
  assert(xfer_size <= VTA_ACC_BUFF_DEPTH);

  // Initialize instruction buffer
  VTAGenericInsn *insn_buf =
      static_cast<VTAGenericInsn *>(allocBuffer(sizeof(VTAGenericInsn) * ins_size));

  // Load uops
  insn_buf[0] = getLoadStoreInsn(VTA_OPCODE_LOAD, VTA_MEM_ID_UOP, uop_size);
  // Load data block
  insn_buf[1] = getLoadStoreInsn(VTA_OPCODE_LOAD, VTA_MEM_ID_ACC, xfer_size);
  // Perform GEMM
  insn_buf[2] = getGEMMInsn(uop_size, true);
  // Store data block
  insn_buf[3] = getLoadStoreInsn(VTA_OPCODE_STORE, VTA_MEM_ID_ACC, xfer_size);

  // Prepare the uop buffer
  VTAGenericUop * uop_buf = getResetUops(batch / VTA_BATCH, out_channels / VTA_BLOCK_OUT);

#if VTA_DEBUG == 1
  printInstruction(ins_size, insn_buf);
  printMicroOp(uop_size, uop_buf);
#endif

  // Initialize input data
  acc_T **inputs = allocInit2dArray<acc_T, VTA_ACC_WIDTH>(batch, out_channels);

  // Reference output
  acc_T **outputs_ref = alloc2dArray<acc_T>(batch, out_channels);
  for (int i = 0; i < batch; i++) {
    for (int j = 0; j < out_channels; j++) {
      outputs_ref[i][j] = 0;
    }
  }

  // Prepare the input buffer
  acc_T *input_buf = static_cast<acc_T *>(allocBuffer(VTA_ACC_ELEM_BYTES * xfer_size));
  pack2dBuffer<acc_T, VTA_ACC_WIDTH>(input_buf,
                                   inputs,
                                   batch,
                                   out_channels,
                                   VTA_BATCH,
                                   VTA_BLOCK_OUT);
  // Prepare the output buffer
  acc_T *output_buf = static_cast<acc_T *>(allocBuffer(VTA_ACC_ELEM_BYTES * xfer_size));

#ifdef NO_SIM
  // Invoke the VTA
  uint64_t t_fpga = vta_alveo(ins_size,
                        insn_buf,
                        uop_buf,
                        NULL,
                        NULL,
                        input_buf,
                        output_buf);
  // Report on timining
  printf("INFO - Synchronization time: %.3lfms\n", static_cast<float>(t_fpga) / 1E6);
#else
  // Invoke the VTA
  vta_alveo(ins_size,
     (volatile insn_T *) insn_buf,
     (volatile uop_T *) uop_buf,
     NULL,
     NULL,
     (volatile acc_vec_T *) input_buf,
     (volatile acc_vec_T *) output_buf);
#endif

  // Unpack output data
  acc_T **outputs = alloc2dArray<acc_T>(batch, out_channels);
  unpack2dBuffer<acc_T, VTA_ACC_WIDTH>(outputs,
                                       output_buf,
                                       batch,
                                       out_channels,
                                       VTA_BATCH,
                                       VTA_BLOCK_OUT);

  // Correctness checks
  int err = 0;
  for (int i = 0; i < batch; i++) {
    for (int j = 0; j < out_channels; j++) {
      if (outputs_ref[i][j] != outputs[i][j]) {
        err++;
#if VTA_DEBUG == 1
        printf("DEBUG - %d, %d: expected 0x%x but got 0x%x\n", i, j,
               static_cast<int>(outputs_ref[i][j]),
               static_cast<int>(outputs[i][j]));
#endif  // VTA_DEBUG
      }
    }
  }

  // Free all allocated arrays
  free2dArray<acc_T>(inputs, batch, out_channels);
  free2dArray<acc_T>(outputs_ref, batch, out_channels);
  free2dArray<acc_T>(outputs, batch, out_channels);
  freeBuffer(insn_buf);
  freeBuffer(uop_buf);
  freeBuffer(input_buf);
  freeBuffer(output_buf);

  if (err == 0) {
    printf("INFO - Reset test successful!\n");
    return 0;
  } else {
    printf("INFO - Reset test failed, got %d errors!\n", err);
    return -1;
  }
}

int fc_test(int batch, int in_channels, int out_channels) {
  // Some assertions
  assert(in_channels % VTA_BLOCK_IN == 0);
  assert(out_channels % VTA_BLOCK_OUT == 0);
  assert(batch % VTA_BATCH == 0);

  printf("=====================================================================================\n");
  printf("INFO - FC test: batch=%d, in_channels=%d, out_channels=%d\n",
         batch, in_channels, out_channels);

  // Derive number of elements that need to be loaded/stored
  int ins_size = 6;
  int uop_size = batch / VTA_BATCH * in_channels / VTA_BLOCK_IN * out_channels / VTA_BLOCK_OUT;
  int inp_size = batch / VTA_BATCH * in_channels / VTA_BLOCK_IN;
  int wgt_size = in_channels / VTA_BLOCK_IN * out_channels / VTA_BLOCK_OUT;
  int out_size = batch / VTA_BATCH * out_channels / VTA_BLOCK_OUT;
  // Make sure we don't exceed buffer bounds
#if VTA_DEBUG == 1
  printf("INFO - uop size = %d/%d\n", uop_size, VTA_UOP_BUFF_DEPTH);
  printf("INFO - input size = %d/%d\n", inp_size, VTA_INP_BUFF_DEPTH);
  printf("INFO - weight size = %d/%d\n", wgt_size, VTA_WGT_BUFF_DEPTH);
  printf("INFO - out size = %d/%d\n", out_size, VTA_ACC_BUFF_DEPTH);
#endif
  assert(uop_size <= VTA_UOP_BUFF_DEPTH);
  assert(inp_size <= VTA_INP_BUFF_DEPTH);
  assert(wgt_size <= VTA_WGT_BUFF_DEPTH);
  assert(out_size <= VTA_ACC_BUFF_DEPTH);

  // Initialize instruction buffer
  VTAGenericInsn *insn_buf =
      static_cast<VTAGenericInsn *>(allocBuffer(sizeof(VTAGenericInsn) * ins_size));

  // Load uops
  insn_buf[0] = getLoadStoreInsn(VTA_OPCODE_LOAD, VTA_MEM_ID_UOP, uop_size);
  // Load bias block
  insn_buf[1] = getLoadStoreInsn(VTA_OPCODE_LOAD, VTA_MEM_ID_ACC, out_size);
  // Load input block
  insn_buf[2] = getLoadStoreInsn(VTA_OPCODE_LOAD, VTA_MEM_ID_INP, inp_size);
  // Load weight block
  insn_buf[3] = getLoadStoreInsn(VTA_OPCODE_LOAD, VTA_MEM_ID_WGT, wgt_size);
  // Perform GEMM
  insn_buf[4] = getGEMMInsn(uop_size, false);
  // Store output block
  insn_buf[5] = getLoadStoreInsn(VTA_OPCODE_STORE, VTA_MEM_ID_ACC, out_size);

  // Prepare the uop buffer
  VTAGenericUop * uop_buf = getGEMMUops(batch / VTA_BATCH,
                                    in_channels / VTA_BLOCK_IN,
                                    out_channels / VTA_BLOCK_OUT);

#if VTA_DEBUG == 1
  printInstruction(ins_size, insn_buf);
  printMicroOp(uop_size, uop_buf);
#endif

  // Initialize inputs
  inp_T **inputs = allocInit2dArray<inp_T, VTA_INP_WIDTH>(batch, in_channels);
  // Initialize weights
  wgt_T **weights = allocInit2dArray<wgt_T, VTA_WGT_WIDTH>(out_channels, in_channels);
  // Initialize biases
  acc_T **biases = allocInit2dArray<acc_T, VTA_ACC_WIDTH>(batch, out_channels);

  // Reference output
  acc_T **outputs_ref = alloc2dArray<acc_T>(batch, out_channels);
  for (int i = 0; i < batch; i++) {
    for (int j = 0; j < out_channels; j++) {
      acc_T sum = biases[i][j];
      for (int k = 0; k < in_channels; k++) {
        sum += (acc_T) (inputs[i][k] * weights[j][k]);
      }
      // Set
      outputs_ref[i][j] = (acc_T) sum;
    }
  }

  // Prepare the input buffer
  inp_T *input_buf = static_cast<inp_T *>(allocBuffer(VTA_INP_ELEM_BYTES * inp_size));
  pack2dBuffer<inp_T, VTA_INP_WIDTH>(input_buf,
                                     inputs,
                                     batch,
                                     in_channels,
                                     VTA_BATCH,
                                     VTA_BLOCK_IN);
  // Prepare the weight buffer
  wgt_T *weight_buf = static_cast<wgt_T *>(allocBuffer(VTA_WGT_ELEM_BYTES * wgt_size));
  pack2dBuffer<wgt_T, VTA_WGT_WIDTH>(weight_buf,
                                     weights,
                                     out_channels,
                                     in_channels,
                                     VTA_BLOCK_OUT,
                                     VTA_BLOCK_IN);
  // Prepare the bias buffer
  acc_T *bias_buf = static_cast<acc_T *>(allocBuffer(VTA_ACC_ELEM_BYTES * out_size));
  pack2dBuffer<acc_T, VTA_ACC_WIDTH>(bias_buf,
                                     biases,
                                     batch,
                                     out_channels,
                                     VTA_BATCH,
                                     VTA_BLOCK_OUT);
  // Prepare the output buffer
  acc_T *output_buf = static_cast<acc_T *>(allocBuffer(VTA_ACC_ELEM_BYTES * out_size));

#ifdef NO_SIM
  // Invoke the VTA
  uint64_t t_fpga = vta_alveo(ins_size,
                        insn_buf,
                        uop_buf,
                        input_buf,
                        weight_buf,
                        bias_buf,
                        output_buf);
  // Report on timining
  printf("INFO - Synchronization time: %.3lfms\n", static_cast<float>(t_fpga) / 1E6);
  printf("INFO - Throughput: %.3lfGOPs/s\n",
         static_cast<float>(uop_size) * VTA_BATCH * VTA_BLOCK_IN * VTA_BLOCK_OUT * 2 / t_fpga);
#else
  // Invoke the VTA
  vta_alveo(ins_size,
     (volatile insn_T *) insn_buf,
     (volatile uop_T *) uop_buf,
     (volatile inp_vec_T *) input_buf,
     (volatile wgt_vec_T *) weight_buf,
     (volatile acc_vec_T *) bias_buf,
     (volatile acc_vec_T *) output_buf);
#endif

  // Unpack output data
  acc_T **outputs = alloc2dArray<acc_T>(batch, out_channels);
  unpack2dBuffer<acc_T, VTA_ACC_WIDTH>(outputs,
                                       output_buf,
                                       batch,
                                       out_channels,
                                       VTA_BATCH,
                                       VTA_BLOCK_OUT);

  // Correctness checks
  int err = 0;
  for (int i = 0; i < batch; i++) {
    for (int j = 0; j < out_channels; j++) {
      if (outputs_ref[i][j] != outputs[i][j]) {
        err++;
#if VTA_DEBUG == 1
        printf("DEBUG - %d, %d: expected 0x%x but got 0x%x\n", i, j,
               static_cast<int>(outputs_ref[i][j]),
               static_cast<int>(outputs[i][j]));
#endif  // VTA_DEBUG
      }
    }
  }

  // Free all allocated arrays
  free2dArray<inp_T>(inputs, batch, in_channels);
  free2dArray<wgt_T>(weights, out_channels, in_channels);
  free2dArray<acc_T>(biases, batch, out_channels);
  free2dArray<acc_T>(outputs_ref, batch, out_channels);
  free2dArray<acc_T>(outputs, batch, out_channels);
  freeBuffer(insn_buf);
  freeBuffer(uop_buf);
  freeBuffer(input_buf);
  freeBuffer(weight_buf);
  freeBuffer(bias_buf);
  freeBuffer(output_buf);

  if (err == 0) {
    printf("INFO - FC test successful!\n");
    return 0;
  } else {
    printf("INFO - FC test failed, got %d errors!\n", err);
    return -1;
  }
}


int conv2d_test(int batch, int height, int width, int kheight, int kwidth,
    int in_channels, int out_channels) {
  // Some assertions
  assert(in_channels % VTA_BLOCK_IN == 0);
  assert(out_channels % VTA_BLOCK_OUT == 0);
  assert(batch % VTA_BATCH == 0);

  printf("=====================================================================================\n");
  printf("INFO - 2D Convolution test: batch=%d, height=%d, width=%d\n", batch, height, width);
  printf("                            kheight=%d, kwidth=%d\n", kheight, kwidth);
  printf("                            in_channels=%d, out_channels=%d\n", in_channels, out_channels);

  // Derive padding
  int ypad = (kheight - 1) / 2;
  int xpad = (kwidth - 1) / 2;
  // Derive number of elements that need to be loaded/stored
  int ins_size = 6;
  int uop_size = batch / VTA_BATCH * in_channels / VTA_BLOCK_IN * out_channels / VTA_BLOCK_OUT *
      (height - kheight + 1) * (width - kwidth + 1) * kheight * kwidth;
  int inp_size = batch / VTA_BATCH * in_channels / VTA_BLOCK_IN * height * width;
  int wgt_size = in_channels / VTA_BLOCK_IN * out_channels / VTA_BLOCK_OUT * kheight * kwidth;
  int out_size = batch / VTA_BATCH * out_channels / VTA_BLOCK_OUT * height * width;
  // Make sure we don't exceed buffer bounds
#if VTA_DEBUG == 1
  printf("INFO - uop size = %d/%d\n", uop_size, VTA_UOP_BUFF_DEPTH);
  printf("INFO - input size = %d/%d\n", inp_size, VTA_INP_BUFF_DEPTH);
  printf("INFO - weight size = %d/%d\n", wgt_size, VTA_WGT_BUFF_DEPTH);
  printf("INFO - out size = %d/%d\n", out_size, VTA_ACC_BUFF_DEPTH);
#endif
  assert(uop_size <= VTA_UOP_BUFF_DEPTH);
  assert(inp_size <= VTA_INP_BUFF_DEPTH);
  assert(wgt_size <= VTA_WGT_BUFF_DEPTH);
  assert(out_size <= VTA_ACC_BUFF_DEPTH);

  // Initialize instruction buffer
  VTAGenericInsn *insn_buf =
      static_cast<VTAGenericInsn *>(allocBuffer(sizeof(VTAGenericInsn) * ins_size));

  // Load uops
  insn_buf[0] = getLoadStoreInsn(VTA_OPCODE_LOAD, VTA_MEM_ID_UOP, uop_size);
  // Load bias block
  insn_buf[1] = getLoadStoreInsn(VTA_OPCODE_LOAD, VTA_MEM_ID_ACC, out_size);
  // Load input block
  insn_buf[2] = getLoadStoreInsn(VTA_OPCODE_LOAD, VTA_MEM_ID_INP, inp_size);
  // Load weight block
  insn_buf[3] = getLoadStoreInsn(VTA_OPCODE_LOAD, VTA_MEM_ID_WGT, wgt_size);
  // Perform GEMM
  insn_buf[4] = getGEMMInsn(uop_size, false);
  // Store output block
  insn_buf[5] = getLoadStoreInsn(VTA_OPCODE_STORE, VTA_MEM_ID_ACC, out_size);

  // Prepare the uop buffer
  VTAGenericUop * uop_buf = getConv2dUops(batch / VTA_BATCH,
                                      height,
                                      width,
                                      kheight,
                                      kwidth,
                                      in_channels / VTA_BLOCK_IN,
                                      out_channels / VTA_BLOCK_OUT);

#if VTA_DEBUG == 1
  printInstruction(ins_size, insn_buf);
  printMicroOp(uop_size, uop_buf);
#endif

  // Initialize inputs
  inp_T ****inputs = allocInit4dArray<inp_T, VTA_INP_WIDTH>(batch, height, width, in_channels);
  // Initialize weights
  wgt_T ****weights = allocInit4dArray<wgt_T, VTA_WGT_WIDTH>(out_channels, kheight, kwidth, in_channels);
  // Initialize biases
  acc_T ****biases = allocInit4dArray<acc_T, VTA_ACC_WIDTH>(batch, height, width, out_channels);

  // Reference output
  acc_T ****outputs_ref = alloc4dArray<acc_T>(batch, height, width, out_channels);
  for (int b = 0; b < batch; b++) {
    for (int y = ypad; y < height - ypad; y++) {
      for (int x = xpad; x < width - xpad; x++) {
        for (int oc = 0; oc < out_channels; oc++) {
          acc_T sum = biases[b][y][x][oc];
          for (int dy = 0 - ypad; dy < ypad + 1; dy++) {
            for (int dx = 0 - xpad; dx < xpad + 1; dx++) {
              for (int ic = 0; ic < in_channels; ic++) {
                sum += (acc_T) (inputs[b][y+dy][x+dx][ic] * weights[oc][dy+ypad][dx+xpad][ic]);
              }
            }
          }
          // Set
          outputs_ref[b][y][x][oc] = sum;
        }
      }
    }
  }

  // Prepare the input buffer
  inp_T *input_buf = static_cast<inp_T *>(allocBuffer(VTA_INP_ELEM_BYTES * inp_size));
  pack4dBuffer<inp_T, VTA_INP_WIDTH>(input_buf,
                                     inputs,
                                     batch,
                                     height,
                                     width,
                                     in_channels,
                                     VTA_BATCH,
                                     VTA_BLOCK_IN);
  // Prepare the weight buffer
  wgt_T *weight_buf = static_cast<wgt_T *>(allocBuffer(VTA_WGT_ELEM_BYTES * wgt_size));
  pack4dBuffer<wgt_T, VTA_WGT_WIDTH>(weight_buf,
                                     weights,
                                     out_channels,
                                     kheight,
                                     kwidth,
                                     in_channels,
                                     VTA_BLOCK_OUT,
                                     VTA_BLOCK_IN);
  // Prepare the bias buffer
  acc_T *bias_buf = static_cast<acc_T *>(allocBuffer(VTA_ACC_ELEM_BYTES * out_size));
  pack4dBuffer<acc_T, VTA_ACC_WIDTH>(bias_buf,
                                     biases,
                                     batch,
                                     height,
                                     width,
                                     out_channels,
                                     VTA_BATCH,
                                     VTA_BLOCK_OUT);
  // Prepare the output buffer
  acc_T *output_buf = static_cast<acc_T *>(allocBuffer(VTA_ACC_ELEM_BYTES * out_size));

#ifdef NO_SIM
  // Invoke the VTA
  uint64_t t_fpga = vta_alveo(ins_size,
                        insn_buf,
                        uop_buf,
                        input_buf,
                        weight_buf,
                        bias_buf,
                        output_buf);
  // Report on timining
  printf("INFO - Synchronization time: %.3lfms\n", static_cast<float>(t_fpga) / 1E6);
  printf("INFO - Throughput: %.3lfGOPs/s\n",
         static_cast<float>(uop_size) * VTA_BATCH * VTA_BLOCK_IN * VTA_BLOCK_OUT * 2 / t_fpga);
#else
  // Invoke the VTA
  vta_alveo(ins_size,
     (volatile insn_T *) insn_buf,
     (volatile uop_T *) uop_buf,
     (volatile inp_vec_T *) input_buf,
     (volatile wgt_vec_T *) weight_buf,
     (volatile acc_vec_T *) bias_buf,
     (volatile acc_vec_T *) output_buf);
#endif

  // Unpack output data
  acc_T ****outputs = alloc4dArray<acc_T>(batch, height, width, out_channels);
  unpack4dBuffer<acc_T, VTA_ACC_WIDTH>(outputs,
                                       output_buf,
                                       batch,
                                       height,
                                       width,
                                       out_channels,
                                       VTA_BATCH,
                                       VTA_BLOCK_OUT);

  // Correctness checks
  // Note we ignore the boundary elements
  int err = 0;
  for (int b = 0; b < batch; b++) {
    for (int y = ypad; y < height - ypad; y++) {
      for (int x = xpad; x < width - xpad; x++) {
        for (int oc = 0; oc < out_channels; oc++) {
          if (outputs_ref[b][y][x][oc] != outputs[b][y][x][oc]) {
            err++;
#if VTA_DEBUG == 1
            printf("DEBUG - %d, %d, %d, %d: expected 0x%x but got 0x%x\n",
                   b, y, x, oc,
                   static_cast<int>(outputs_ref[b][y][x][oc]),
                   static_cast<int>(outputs[b][y][x][oc]));
#endif  // VTA_DEBUG == 1
          }
        }
      }
    }
  }

  // Free all allocated arrays
  free4dArray<inp_T>(inputs, batch, height, width, in_channels);
  free4dArray<wgt_T>(weights, out_channels, kheight, kwidth, in_channels);
  free4dArray<acc_T>(biases, batch, height, width, out_channels);
  free4dArray<acc_T>(outputs_ref, batch, height, width, out_channels);
  free4dArray<acc_T>(outputs, batch, height, width, out_channels);
  freeBuffer(insn_buf);
  freeBuffer(uop_buf);
  freeBuffer(input_buf);
  freeBuffer(weight_buf);
  freeBuffer(bias_buf);
  freeBuffer(output_buf);

  if (err == 0) {
    printf("INFO - 2D Convolution test successful!\n");
    return 0;
  } else {
    printf("INFO - 2D Convolution test failed, got %d errors!\n", err);
    return -1;
  }
}
