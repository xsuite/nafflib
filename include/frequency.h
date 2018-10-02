#ifndef __NAFFLIB_FREQUENCY_H__
#define __NAFFLIB_FREQUENCY_H__

#include <complex.h>
#include <math.h>
#include "signal_processing.h"
#include "fft.h"
#include "brent.h"
#include "windows.h"

struct merit_args
{
    size_t N;
    double _Complex* window;
    double _Complex* signal;

};

typedef struct merit_args merit_args;

double minus_magnitude_fourier_integral(double frequency, const merit_args* S);
double get_f1(double _Complex* signal, size_t N, double order, double fft_estimate);
void py_f1(double _Complex* signal, size_t N, double order, double fft_estimate, double* tune);
void pyget_f1(double _Complex* signal, size_t N, double order, double *tune);
void get_tune(double _Complex* signal, size_t N, double order, double *tune);
#endif
