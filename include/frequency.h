#ifndef __NAFFLIB_FREQUENCY_H__
#define __NAFFLIB_FREQUENCY_H__

#include <complex.h>
#include <math.h>
#include "signal_processing.h"
#include "fft.h"
#include "brent.h"
#include "windows.h"

//#define MAX_FREQS 1000
//
struct merit_args
{
    size_t N;
    double _Complex* window;
    double _Complex* signal;
    
//    double _Complex amplitude[MAX_FREQS];
//    double  frequency[MAX_FREQS];
//    int frequency_counter;
//
};

typedef struct merit_args merit_args;

double minus_magnitude_fourier_integral(double frequency, const merit_args* S);
double get_f1(double _Complex* signal, size_t N, double order, double fft_estimate);
size_t interpolating_size(size_t N);
void use_interpolating_integral(size_t N, double _Complex* window);
void get_f_neg(double _Complex* signal, size_t N, double order, double* frequencies, double _Complex* amplitudes, double _Complex* negative_amplitudes, size_t n_freqs);
void py_f1(double _Complex* signal, size_t N, double order, double fft_estimate, double* tune);
void pyget_f1(double _Complex* signal, size_t N, double order, double *tune);
void get_tune(double _Complex* signal, size_t N, double order, double *tune);
#endif
