#ifndef _NAFFLIB_BRENT_H__
#define _NAFFLIB_BRENT_H__

#include "stdio.h"
#include "math.h"
#include "frequency.h"

double brent_minimize(double (*f)(double, const merit_args*), double min, double max,const merit_args* S);
#endif
