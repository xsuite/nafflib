#ifndef __NAFFLIB_SIGNAL_PROCESSING_H__
#define __NAFFLIB_SIGNAL_PROCESSING_H__

#include <complex.h>
#include <stdio.h>
#include <math.h>

#define pi 3.141592653589793238462643383279

#ifdef COMPILE_WITH_CRLIBM
#include "crlibm.h"
#define cos cos_rn
#define sin sin_rn
#endif

double _Complex inner_product(const double _Complex* signal, double amplitude, double frequency,const double _Complex* window, size_t N);
void strip_DC(double _Complex* signal, size_t N);

double _Complex frequency_project(double f0, double _Complex *amps, double *freqs, size_t n_frequencies, double _Complex *window, size_t N);
void subtract_frequency(double _Complex* signal, double _Complex** normal_amps, double* frequencies, double _Complex* amplitudes, size_t i, double _Complex* window, size_t N);
void remove_component( double _Complex* signal, double _Complex* amps, double* frequencies, size_t n_freqs, size_t N);
double _Complex signal_project(double _Complex* signal, double _Complex* amps, double* freqs, size_t n_frequencies, double _Complex* window, size_t N);

#endif
