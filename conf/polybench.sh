#!/bin/bash

export APPDIR="./polybench"
export BINARYS="datamining/correlation/correlation|LARGE
datamining/covariance/covariance|LARGE
linear-algebra/blas/gemm/gemm|LARGE
linear-algebra/blas/gemver/gemver|LARGE
linear-algebra/blas/gesummv/gesummv|LARGE
linear-algebra/blas/symm/symm|LARGE
linear-algebra/blas/syr2k/syr2k|LARGE
linear-algebra/blas/syrk/syrk|LARGE
linear-algebra/blas/trmm/trmm|LARGE
linear-algebra/kernels/2mm/2mm|LARGE
linear-algebra/kernels/3mm/3mm|LARGE
linear-algebra/kernels/atax/atax|LARGE
linear-algebra/kernels/bicg/bicg|LARGE
linear-algebra/kernels/doitgen/doitgen|LARGE
linear-algebra/kernels/mvt/mvt|LARGE
linear-algebra/solvers/cholesky/cholesky|LARGE
linear-algebra/solvers/durbin/durbin|LARGE
linear-algebra/solvers/gramschmidt/gramschmidt|LARGE
linear-algebra/solvers/lu/lu|LARGE
linear-algebra/solvers/ludcmp/ludcmp|LARGE
linear-algebra/solvers/trisolv/trisolv|LARGE
medley/deriche/deriche|LARGE
medley/floyd-warshall/floyd-warshall|MEDIUM
medley/nussinov/nussinov|LARGE
stencils/adi/adi|LARGE
stencils/fdtd-2d/fdtd-2d|LARGE
stencils/heat-3d/heat-3d|LARGE
stencils/jacobi-1d/jacobi-1d|LARGE
stencils/jacobi-2d/jacobi-2d|LARGE
stencils/seidel-2d/seidel-2d|MEDIUM"
export INPUT=""
export NumRunsTEST=1
export NumRunsBEST=10
export MAXTIME="2m"
export RUNSDE="yes"
export RUNPCM="no"
export RUNVTUNE="no"
export TESTCONF="1"
export BESTCONF="1"
