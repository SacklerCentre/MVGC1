Compiling MVGC C source
-----------------------

There is a gcc-compatible Makefile in this directory. On linux x86_64 with gcc, the mex command will look something like:

    mex -O -largeArrayDims CFLAGS='$CFLAGS -Wall -Werror -std=c99 -O3' -lut -lm genvar_mex.c -o ../mex/genvar_mex.mexa64

On other platforms your starting point is

    mex genvar_mex.c -o ../mex/genvar_mex.<arch>
