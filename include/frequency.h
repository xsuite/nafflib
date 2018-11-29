#ifndef __NAFFLIB_FREQUENCY_H__
#define __NAFFLIB_FREQUENCY_H__

#include <complex.h>
#include <math.h>
#include "signal_processing.h"
#include "naffargs.h"
#include "fft.h"
#include "brent.h"
#include "windows.h"

double minus_magnitude_fourier_integral(double frequency, const merit_args* S);
double get_f1(double _Complex* signal, size_t N, double order, double fft_estimate);
size_t interpolating_size(size_t N);
void use_interpolating_integral_hardy(size_t N, double _Complex* window);
void use_interpolating_integral(size_t N, double _Complex* window);

double get_q(double _Complex* signal, size_t N, double order, int interpolate_integral);
void get_f_all(double _Complex* signal, size_t N, double order,  double* frequencies, double _Complex* amplitudes, size_t n_freqs);
void get_f_neg(double _Complex* signal, size_t N, double order, double* frequencies, double _Complex* amplitudes, double _Complex* negative_amplitudes, size_t n_freqs);

#endif
