#Number of compute units 
NK:= 2
#Number of parallel points
PP:= 96

CXXFLAGS += -D __USE_OPENCL__ -DNUM_CU=$(NK)
#CXXFLAGS += -D __USE_OPENCL__ -DNUM_CU=$(NK) -D VERIFY_USING_CMODEL
CLFLAGS += -DPARALLEL_POINTS=$(PP)

