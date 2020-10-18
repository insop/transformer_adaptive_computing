#include <algorithm>
#include <iostream>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#include "xcl2.hpp"

// HBM Banks requirements
#define MAX_HBM_BANKCOUNT 32
#define BANK_NAME(n) n | XCL_MEM_TOPOLOGY
const int bank[MAX_HBM_BANKCOUNT] = {
    BANK_NAME(0),  BANK_NAME(1),  BANK_NAME(2),  BANK_NAME(3),  BANK_NAME(4),
    BANK_NAME(5),  BANK_NAME(6),  BANK_NAME(7),  BANK_NAME(8),  BANK_NAME(9),
    BANK_NAME(10), BANK_NAME(11), BANK_NAME(12), BANK_NAME(13), BANK_NAME(14),
    BANK_NAME(15), BANK_NAME(16), BANK_NAME(17), BANK_NAME(18), BANK_NAME(19),
    BANK_NAME(20), BANK_NAME(21), BANK_NAME(22), BANK_NAME(23), BANK_NAME(24),
    BANK_NAME(25), BANK_NAME(26), BANK_NAME(27), BANK_NAME(28), BANK_NAME(29),
    BANK_NAME(30), BANK_NAME(31)};

int main(int argc, char *argv[]) {
  if (argc != 2) {
    printf("Usage: %s <XCLBIN> \n", argv[0]);
    return -1;
  }

  if (xcl::is_emulation()) {
  }

  /* Allocate space in each one of 16 HBM banks */
  /* Load (3072,1024) weights into HBM in Device Global memory*/

  /* Load (1024,14) vector from Host Local into separate HBM in Device Global memory */
  /* Run both control_1() and wb_1() kernels */
  /* Load (1024,14) results vector from Device Global memory back to Host Local */


  std::string binaryFile = argv[1];
  cl_int err;
  cl::CommandQueue q[2];
  cl::Kernel krnl_control;
  cl::Kernel krnl_wb;
  cl::Context context;

#if 0
  std::vector<float, aligned_allocator<float>> source_w[16];
  std::vector<float, aligned_allocator<float>> source_v(sizeof(signed short)*14*1024);
  std::vector<float, aligned_allocator<float>> source_hw_wb_results(sizeof(signed short)*14*1024);

  for(int i = 0; i < 16; i++) {
    source_w[i].resize(sizeof(signed short)*3*1024*1024/16);
  }

  // Create the test data
  for(int i = 0; i < 16; i++) {
    std::generate(source_w[i].begin(), source_w[i].end(), std::rand);
  }
  std::generate(source_v.begin(), source_v.end(), std::rand);

  // Initializing output vectors to zero
  std::fill(source_hw_wb_results.begin(), source_hw_wb_results.end(), 0);
#endif

  // OPENCL HOST CODE AREA START
  // The get_xil_devices will return vector of Xilinx Devices
  auto devices = xcl::get_xil_devices();

  // read_binary_file() command will find the OpenCL binary file created using
  // the  V++ compiler load into OpenCL Binary and return pointer to file buffer.
  auto fileBuf = xcl::read_binary_file(binaryFile);
  cl::Program::Binaries bins{{fileBuf.data(), fileBuf.size()}};
  bool valid_device = false;
  for (unsigned int i = 0; i < devices.size(); i++) {
    auto device = devices[i];
    // Creating Context and Command Queue for selected Device
    OCL_CHECK(err, context = cl::Context(device, NULL, NULL, NULL, &err));
    for(int i = 0; i < 2; i++) {
      OCL_CHECK(err, q[i] = cl::CommandQueue(context, device,
                                        CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE |
                                            CL_QUEUE_PROFILING_ENABLE,
                                        &err));
    }
    std::cout << "Trying to program device[" << i
              << "]: " << device.getInfo<CL_DEVICE_NAME>() << std::endl;
    cl::Program program(context, {device}, bins, NULL, &err);
    if (err != CL_SUCCESS) {
      std::cout << "Failed to program device[" << i << "] with xclbin file!\n";
    } else {
      std::cout << "Device[" << i << "]: program successful!\n";
      // Creating Kernel object using Compute unit names

      std::string krnl_name_control = "vta_alveo";
      printf("Creating a kernel [%s] for CU\n", krnl_name_control.c_str());
      OCL_CHECK(err, krnl_control = cl::Kernel(program, krnl_name_control.c_str(), &err));

#if 0
      std::string krnl_name_wb = "wb:{wb_1}";
      printf("Creating a kernel [%s] for CU\n", krnl_name_wb.c_str());
      OCL_CHECK(err, krnl_wb = cl::Kernel(program, krnl_name_wb.c_str(), &err));
#endif
      valid_device = true;
      break; // we break because we found a valid device
    }
  }
  if (!valid_device) {
    std::cout << "Failed to program any device found, exit!\n";
    exit(EXIT_FAILURE);
  }

#if 0
  std::vector<cl_mem_ext_ptr_t> inBufExtw(16);
  std::vector<cl_mem_ext_ptr_t> inBufExtv(1);
  std::vector<cl_mem_ext_ptr_t> outBufExtwb(1);

  std::vector<cl::Buffer> buffer_inputw(16);
  std::vector<cl::Buffer> buffer_inputv(1);
  std::vector<cl::Buffer> buffer_output_wb(1);

  // For Allocating Buffer to specific Global Memory Bank, user has to use
  // cl_mem_ext_ptr_t
  // and provide the Banks
  for(int i = 0; i < 16; i++) {
    inBufExtw[i].obj = source_w[i].data();
    inBufExtw[i].param = 0;
    inBufExtw[i].flags = bank[i * 2];
  }

  inBufExtv[0].obj = source_v.data();
  inBufExtv[0].param = 0;
  inBufExtv[0].flags = bank[31];

  outBufExtwb[0].obj = source_hw_wb_results.data();
  outBufExtwb[0].param = 0;
  outBufExtwb[0].flags = bank[31];

  // These commands will allocate memory on the FPGA and copy the data across
  // The cl::Buffer objects can be used to reference the memory locations on the device.
  // The Weights are spread across 16 even-numbered HBM banks and vector lives in bank 15.
  for(int i = 0; i < 16; i++) {
    OCL_CHECK(err, buffer_inputw[i] = cl::Buffer(
                      context, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX | CL_MEM_USE_HOST_PTR,
                      sizeof(signed short) * 3 * 1024 * 1024/16, &inBufExtw[i], &err));
    OCL_CHECK(err, err = q[0].enqueueMigrateMemObjects({buffer_inputw[i]}, 0 /* 0 means from host*/));
  }
  OCL_CHECK(err, buffer_inputv[0] = cl::Buffer(
                      context, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX | CL_MEM_USE_HOST_PTR,
                      sizeof(signed short) * 1024 * 14, &inBufExtv[0], &err));
  OCL_CHECK(err, err = q[0].enqueueMigrateMemObjects({buffer_inputv[0]}, 0 /* 0 means from host*/));

  OCL_CHECK(err, buffer_output_wb[0] = cl::Buffer(
                      context, CL_MEM_WRITE_ONLY | CL_MEM_EXT_PTR_XILINX | CL_MEM_USE_HOST_PTR,
                      sizeof(signed short) * 1024 * 14, &outBufExtwb[0], &err));
  OCL_CHECK(err, err = q[0].enqueueMigrateMemObjects({buffer_output_wb[0]}, 0 /* 0 means from host*/));
  // Copy input data to Device Global Memory
  q[0].finish();

  int outshiftscale = 0;
  OCL_CHECK(err, err = krnl_wb.setArg(0, buffer_output_wb[0]));
  OCL_CHECK(err, err = krnl_wb.setArg(1, outshiftscale));

  double kernel_time_in_sec = 0;
  std::chrono::duration<double> kernel_time(0);
  auto kernel_start = std::chrono::high_resolution_clock::now();
  // Setting the control Arguments
  for(int i = 0; i < 16; i++) {
    OCL_CHECK(err, err = krnl_control.setArg(i, buffer_inputw[i]));
  }
  OCL_CHECK(err, err = krnl_control.setArg(16, buffer_inputv[0]));

  
  // Invoking the kernel
  OCL_CHECK(err, err = q[0].enqueueTask(krnl_control));
  OCL_CHECK(err, err = q[1].enqueueTask(krnl_wb));
  q[0].finish();
  q[1].finish();

  auto kernel_end = std::chrono::high_resolution_clock::now();
  kernel_time = std::chrono::duration<double>(kernel_end - kernel_start);
  kernel_time_in_sec = kernel_time.count();

  // Copy Result from Device Global Memory to Host Local Memory
  OCL_CHECK(err, err = q[0].enqueueMigrateMemObjects({buffer_output_wb[0]}, CL_MIGRATE_MEM_OBJECT_HOST));
  q[0].finish();

  std::cout << "kernel time = " << (kernel_time_in_sec*1000000) << " us" << std::endl;
  // OPENCL HOST CODE AREA ENDS

  int match = 1;
  /* Need to do some verification of the result*/
  std::cout << (match ? "TEST PASSED" : "TEST FAILED") << std::endl;
  return (match ? EXIT_SUCCESS : EXIT_FAILURE);
#endif
  return 0;
}
