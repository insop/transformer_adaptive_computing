#include <algorithm>
#include <iostream>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#include "xcl2.hpp"

#include "test_lib.h"

static struct ocl_ctx OCL_CTX;
struct ocl_ctx *OCL_CTX_p = &OCL_CTX;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    printf("Usage: %s <XCLBIN> \n", argv[0]);
    return -1;
  }

  if (xcl::is_emulation()) {
  }

  std::string binaryFile = argv[1];
  cl_int err;
  //cl::CommandQueue q;
  //cl::Kernel krnl_vta_alveo;
  //cl::Context context;

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
    OCL_CHECK(err, OCL_CTX_p->context = cl::Context(device, NULL, NULL, NULL, &err));
    OCL_CHECK(err, OCL_CTX_p->q = cl::CommandQueue(OCL_CTX_p->context, device,
                                      //CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE |
                                          CL_QUEUE_PROFILING_ENABLE,
                                        &err));
    std::cout << "Trying to program device[" << i
              << "]: " << device.getInfo<CL_DEVICE_NAME>() << std::endl;

    cl::Program program(OCL_CTX_p->context, {device}, bins, NULL, &err);
    if (err != CL_SUCCESS) {
      std::cout << "Failed to program device[" << i << "] with xclbin file!\n";
    } else {
      std::cout << "Device[" << i << "]: program successful!\n";
      // Creating Kernel object using Compute unit names

      std::string krnl_name_control = "vta_alveo";
      printf("Creating a kernel [%s] for CU\n", krnl_name_control.c_str());
      OCL_CHECK(err, OCL_CTX_p->kernel = cl::Kernel(program, krnl_name_control.c_str(), &err));

      valid_device = true;
      break; // we break because we found a valid device
    }
  }
  if (!valid_device) {
    std::cout << "Failed to program any device found, exit!\n";
    exit(EXIT_FAILURE);
  }

#if 1
  printf("%s:%d context 0x%p\n", __func__, __LINE__, OCL_CTX_p->context);
  //OCL_CTX_p->context = context;
  //OCL_CTX_p->q = q;
  //OCL_CTX_p->krnl_vta_alveo = krnl_vta_alveo;

  // Run load/store test
  //int status = mem_test(128, 128);
  int status = mem_test(128, 128);
  assert(status == 0);

  // Run reset test
  status |= reset_test(128, 128);
  assert(status == 0);
  // Run fully connected layer test
  status |= fc_test(64, 128, 128);
  assert(status == 0);

  return status;

#endif

  return 0;
}
