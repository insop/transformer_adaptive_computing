/*!
 *  Copyright (c) 2018 by Contributors
 * \file test_lib.cpp
 * \brief Test library for the VTA design simulation and driver tests.
 */

#ifndef TESTS_HARDWARE_COMMON_TEST_LIB_H_
#define TESTS_HARDWARE_COMMON_TEST_LIB_H_

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>

#include "xcl2.hpp"


// opencl context
struct ocl_ctx {
  cl::Context context;
  cl::CommandQueue q;
  cl::Kernel kernel;
};

extern struct ocl_ctx *OCL_CTX_p;

#define NO_SIM 1

#ifdef NO_SIM

  // #include "../driver/pynq_driver.h"
  #include "./hw_spec.h"

  typedef int8_t wgt_T;
  typedef int8_t inp_T;
  typedef int32_t acc_T;

#if 0
  //uint64_t vta_alveo(uint32_t insn_count,
               //VTAGenericInsn *insns,
               //VTAGenericUop *uops,
               //inp_T *inputs,
               //wgt_T *weights,
               //acc_T *biases,
               //acc_T *outputs);
//
uint64_t vta_alveo(uint32_t insn_count,
             //VTAGenericInsn *insns,
             std::vector<VTAGenericInsn, aligned_allocator<VTAGenericInsn>>  &insns,
             //VTAGenericUop *uops,
             std::vector<VTAGenericUop, aligned_allocator<VTAGenericUop>>  &uops,
             inp_T *inputs,
             wgt_T *weights,
             //acc_T *biases,
             std::vector<acc_T, aligned_allocator<acc_T>> &biases,
             //acc_T *outputs
             std::vector<acc_T, aligned_allocator<acc_T>> &outputs
						);
#endif // 0

#else  // NO_SIM

  #include "./vta_alveo.h"

#endif  // NO_SIM

/*!
* \brief Performs data packing and tiling of 2D array into 1D buffer.
* \param dst Pointer to the packed, and tiled destination 1D array (flattened).
* \param src Pointer to the unpacked source 2D array.
* \param y_size Number of rows.
* \param x_size Number of columns.
* \param y_block Inner tiling along row dimension.
* \param x_block Inner tiling along column dimension.
*/
template <typename T, int T_WIDTH>
void pack2dBuffer(T *dst, T **src, int y_size, int x_size, int y_block, int x_block);

/*!
* \brief Performs data unpacking from a flat 1D buffer into a 2D array.
* \param dst Pointer to the unpacked destination 2D array.
* \param src Pointer to the packed, and tiled source 1D array (flattened).
* \param y_size Number of rows.
* \param x_size Number of columns.
* \param y_block Inner tiling along row dimension.
* \param x_block Inner tiling along column dimension.
*/
template <typename T, int T_WIDTH>
void unpack2dBuffer(T **dst, T *src, int y_size, int x_size, int y_block, int x_block);

/*!
* \brief Performs data packing and tiling of 4D array into 1D buffer.
* \param dst Pointer to the packed, and tiled destination 1D array (flattened).
* \param src Pointer to the unpacked source 4D array.
* \param ax0 Axis 0 length.
* \param ax1 Axis 1 length.
* \param ax2 Axis 2 length.
* \param ax3 Axis 3 length.
* \param ax0_block Outer tiling along axis 0.
* \param ax1_block Inner tiling along axis 3.
*/
template <typename T, int T_WIDTH>
void pack4dBuffer(T *dst, T ****src, int ax0, int ax1, int ax2, int ax3,
    int ax0_block, int ax3_block);

/*!
* \brief Performs data unpacking from a flat 1D buffer into a 4D array.
* \param dst Pointer to the unpacked destination 4D array.
* \param src Pointer to the packed, and tiled source 1D array (flattened).
* \param ax0 Axis 0 length.
* \param ax1 Axis 1 length.
* \param ax2 Axis 2 length.
* \param ax3 Axis 3 length.
* \param ax0_block Outer tiling along axis 0.
* \param ax1_block Inner tiling along axis 3.
*/
template <typename T, int T_WIDTH>
void unpack4dBuffer(T ****dst, T *src, int ax0, int ax1, int ax2, int ax3,
    int ax0_block, int ax3_block);

/*!
* \brief Allocates and initializes a 2D array in the heap.
* \param rows Number of rows.
* \param cols Number of columns.
* \return Pointer to the 2D array.
*/
template <typename T, int T_WIDTH>
T ** allocInit2dArray(int rows, int cols);

/*!
* \brief Allocates a 2D array in the heap.
* \param rows Number of rows.
* \param cols Number of columns.
* \return Pointer to the 2D array.
*/
template <typename T>
T ** alloc2dArray(int rows, int cols);

/*!
* \brief Allocates and initializes a 4D array in the heap.
* \param ax0 Axis 0 length.
* \param ax1 Axis 1 length.
* \param ax2 Axis 2 length.
* \param ax3 Axis 3 length.
* \return Pointer to the 2D array.
*/
template <typename T, int T_WIDTH>
T **** allocInit4dArray(int ax0, int ax1, int ax2, int ax3);

/*!
* \brief Allocates a 4D array in the heap.
* \param ax0 Axis 0 length.
* \param ax1 Axis 1 length.
* \param ax2 Axis 2 length.
* \param ax3 Axis 3 length.
* \return Pointer to the 2D array.
*/
template <typename T>
T **** alloc4dArray(int ax0, int ax1, int ax2, int ax3);

/*!
* \brief Frees a 4D array.
* \param array Pointer to the 4D array to be freed.
* \param ax0 Axis 0 length.
* \param ax1 Axis 1 length.
* \param ax2 Axis 2 length.
* \param ax3 Axis 3 length.
*/
template <typename T>
void free4dArray(T ****array, int ax0, int ax1, int ax2, int ax3);

/*!
* \brief Performs memory allocation in a physically contiguous region of memory.
* \param num_bytes Size of the buffer in bytes.
* \return Pointer to the allocated buffer.
*/
void * allocBuffer(size_t num_bytes);

/*!
* \brief Frees buffer allocated in a physically contiguous region of memory.
* \param buffer Pointer to the buffer to free.
*/
void freeBuffer(void * buffer);

/*!
* \brief Returns a VTA 1D load or store instruction.
* \param opcode Type of operation.
* \param type On-chip memory target.
* \param size Number of elements to load/store.
* \return A VTAGenericInsn for a 1D load or store op.
*/
VTAGenericInsn getLoadStoreInsn(int opcode, int type, int size);

/*!
* \brief Returns a VTA GEMM instruction (can be used for matrix multiply,
*        2D convolution or accum memory reset).
* \param uop_size Size of micro op program in instructions.
* \param reset Reset accumulator memory.
* \return A VTAGenericInsn for a GEMM op.
*/
VTAGenericInsn getGEMMInsn(int uop_size, bool reset);

/*!
* \brief Returns an allocated buffer of VTA micro-ops to implement an accum memory
*        reset operation over a (a, c) matrix.
* \param batch Batch size (a).
* \param out_feat Output features (c).
* \return A VTAGenericUop pointer to an allocated micro-op kernel buffer.
*/
VTAGenericUop * getResetUops(int batch, int out_feat);

/*!
* \brief Returns an allocated buffer of VTA micro-ops to implement a matrix multiplication
*   of size (a, b) x (b, c).
* \param batch Batch size (a).
* \param in_feat Input features (b).
* \param out_feat Output features (c).
* \return A VTAGenericUop pointer to an allocated micro-op kernel buffer.
*/
VTAGenericUop * getGEMMUops(int batch, int in_feat, int out_feat);

/*!
* \brief Returns an allocated buffer of VTA micro-ops to implement a 2d convolution
*        of an input feature map of size (b, h, w, ic), and a kernel map of size
*        (oc, kh, kw, ic), to produce an output feature map of size (b, h, w, oc).
* \param batch Batch size (b).
* \param height Feature map spatial height (h).
* \param width Feature map spatial width (w).
* \param kheight Kernel spatial height (kh).
* \param kwidth Kernel spatial width (kw).
* \param in_feat Input channel depth (ic).
* \param out_feat Output channel depth (oc).
* \return A VTAGenericUop pointer to an allocated micro-op kernel buffer.
*/
VTAGenericUop * getConv2dUops(int batch, int height, int width, int kheight, int kwidth,
    int in_channels, int out_channels);

/*!
* \brief Print out parameters of the VTA design (for debugging purposes).
*/
void printParameters();

/*!
* \brief Print out instruction information (for debugging purposes).
* \param num_insn Number of instructions.
* \param insns Pointer to the instruction buffer.
*/
void printInstruction(int num_insn, VTAGenericInsn *insns);

/*!
* \brief Print out micro-op information (for debugging purposes).
* \param num_insn Number of micro-ops.
* \param insns Pointer to the micro-op buffer.
*/
void printMicroOp(int num_uop, VTAGenericUop *uops);

/*!
* \brief VTA load/store unit test.
*/
int mem_test(int batch, int out_channels);

/*!
* \brief VTA reset accumulator memory unit test.
*/
int reset_test(int batch, int out_channels);

/*!
* \brief VTA FC unit test.
*/
int fc_test(int batch, int in_channels, int out_channels);

/*!
* \brief VTA 2d Convolution unit test.
*/
int conv2d_test(int batch, int height, int width, int kheight, int kwidth,
    int in_channels, int out_channels);

#endif  //  TESTS_HARDWARE_COMMON_TEST_LIB_H_
