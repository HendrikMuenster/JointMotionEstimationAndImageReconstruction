clear all;close all;clc;


mex -largeArrayDims -outdir algorithms/mex -v COMPFLAGS="$COMPFLAGS -openmp" OPTIMFLAGS="$OPTIMFLAGS"  LINKFALGS="$LINKFALGS -openmp" CFLAGS="$CFLAGS -fopenmp" CXXFLAGS="$CXXFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp" 'algorithms/mex/L1TVOpticalFlowCPP.cpp' 'algorithms/mex/tools.cpp'

