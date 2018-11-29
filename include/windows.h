#ifndef __NAFFLIB_WINDOWS_H__
#define __NAFFLIB_WINDOWS_H__


#include <stdlib.h>
#include <math.h>
#include <complex.h> 
#include <assert.h>
#include <stdio.h>

#define pi 3.141592653589793238462643383279

#ifdef COMPILE_WITH_CRLIBM
#include "crlibm.h"
#define log log_rn
#define exp exp_rn
#endif



double cheb_poly(int n, double x);
void cheb_window(double _Complex window[], const size_t N, const double a);
void hann_harm_window(double _Complex window[], const size_t N, const double n);
double calculateFm(size_t m, double sp2, double a, size_t nBar);
void taylorWindow( double _Complex w[], const size_t N, const double NBAR, const double SLL);
void taylor_window(double _Complex window[], const size_t N, const double SLL);
void no_window( double _Complex window[], const size_t N );

#endif
