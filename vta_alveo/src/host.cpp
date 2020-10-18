#include <algorithm>
#include <iostream>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#include "xcl2.hpp"

#include "test_lib.h"

int main(int argc, char *argv[]) {
  if (argc != 2) {
    printf("Usage: %s <XCLBIN> \n", argv[0]);
    return -1;
  }
#if 1

  // Run load/store test
  int status = mem_test(128, 128);

  // Run reset test
  status |= reset_test(128, 128);
  // Run fully connected layer test
  status |= fc_test(64, 128, 128);

  return 0;

#endif

  if (xcl::is_emulation()) {
  }

  std::string binaryFile = argv[1];
  cl_int err;
  cl::CommandQueue q;
  cl::Kernel krnl_control;
  cl::Kernel krnl_wb;
  cl::Context context;

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
    OCL_CHECK(err, q = cl::CommandQueue(context, device,
                                      CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE |
                                          CL_QUEUE_PROFILING_ENABLE,
                                        &err));
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

      valid_device = true;
      break; // we break because we found a valid device
    }
  }
  if (!valid_device) {
    std::cout << "Failed to program any device found, exit!\n";
    exit(EXIT_FAILURE);
  }

  return 0;
}
