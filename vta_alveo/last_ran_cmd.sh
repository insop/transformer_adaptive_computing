
# to bulid hw_emu
make all TARGET=hw_emu DEVICE=xilinx_u50_gen3x16_xdma_201920_3 HOST_ARCH=x86

# to run sw_emu
XCL_EMULATION_MODE=sw_emu ./kmeans -x ./build_dir.sw_emu.xilinx_u50_gen3x16_xdma_201920_3/krnl_kmeans.xclbin -i ./data/100 -c ./data/100.gold_c10 -n 1

# to bulid sw_emu
make all TARGET=sw_emu DEVICE=xilinx_u50_gen3x16_xdma_201920_3 HOST_ARCH=x86
