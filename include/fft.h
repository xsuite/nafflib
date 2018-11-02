#ifndef __NAFFLIB_FFT_H__
#define __NAFFLIB_FFT_H__

#include <stdio.h>
#include <complex.h>
#include <stdlib.h>

#define pi 3.141592653589793238462643383279

void separate(double _Complex* A, size_t N);
void fft2(double _Complex* X, size_t N);
double max_fft_frequency(double _Complex* signal, size_t N);
double _Complex* init_fft(double _Complex* signal, size_t N, size_t *M);

#ifndef COMPILE_WITHOUT_FFTW
#include <fftw3.h>

double max_fftw_frequency_all(double _Complex* signal, size_t N);
double max_fftw_frequency(double _Complex* signal, size_t N);
#endif

#endif
