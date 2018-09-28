#include "signal_processing.h"

double _Complex inner_product(const double _Complex* signal, double amplitude, double frequency,const double _Complex* window, size_t N)
{
    double omega = 2*pi*frequency;
    double _Complex result = 0;
    for( size_t i = N; i--; )
        result += signal[i]*cexp(-I*omega*i)*window[i];
    return amplitude*result/N;
}


