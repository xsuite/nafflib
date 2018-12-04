#ifndef __NAFFLIB_FFT_H__
#define __NAFFLIB_FFT_H__

#include <stdio.h>
#include <complex.h>
#include <stdlib.h>
#include <math.h>

#define pi 3.141592653589793238462643383279

#ifdef COMPILE_WITH_CRLIBM
#include "crlibm.h"
#define sin sin_rn
#endif

size_t reverseBits(size_t num, size_t bitsize);
size_t bitSize(size_t N);
void fft2(double _Complex* data, size_t nn);
double max_fft_frequency(double _Complex* signal, size_t N);
double max_fft_frequency_all(double _Complex* signal, size_t N);
double _Complex* init_fft(double _Complex* signal, size_t N, size_t *M);

#endif
