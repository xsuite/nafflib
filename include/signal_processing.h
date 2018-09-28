#ifndef __NAFFLIB_SIGNAL_PROCESSING_H__
#define __NAFFLIB_SIGNAL_PROCESSING_H__

#include <complex.h>
#include <stdio.h>

#define pi 3.141592653589793238462643383279

double _Complex inner_product(const double _Complex* signal, double amplitude, double frequency,const double _Complex* window, size_t N);
void strip_DC(double _Complex* signal, size_t N);
#endif
