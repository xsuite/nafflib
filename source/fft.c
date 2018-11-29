#include "fft.h"

#include <math.h>

size_t reverseBits(size_t num, size_t bitsize)
{
    size_t reverse_num = 0;
    size_t i = 0;
    for( i = 0; i < bitsize; i++)
    {
        if ((num & (1 << i)))
            reverse_num |= 1 << ((bitsize - 1) - i);
    }
    return reverse_num;
}

size_t bitSize(size_t N)
{
    size_t b = 0;
    while(N != 1)
    {
        N/=2;
        b++;
    }
    return b;
}

void fft2(double _Complex* data, size_t nn)
{
    ////Prepare fft by rearranging data
    //for (size_t i = 0; i < nn; i++)
    //{
    //    size_t j = reverseBits(i,bitSize(nn));
    //    if( j > i)
    //    {
    //        double _Complex temp = data[j];
    //        data[j] = data[i];
    //        data[i] = temp;
    //    }
    //}

    size_t mmax = 2;
    size_t n = nn << 1;
    while( n > mmax)
    {
        size_t mmax2 = mmax >> 1;
        double theta = -2.*pi/mmax;
        double st = sin(0.5*theta);
        double _Complex wp = -2.*st*st + I*sin(theta);
        double _Complex w = 1.;
        for( size_t m = 0; m < mmax2; m++)
        {
            for( size_t i = m; i < nn; i+=mmax)
            {
                size_t j = i+mmax2;
                double _Complex temp = w*data[j];
                data[j] = data[i] - temp;
                data[i] = data[i] + temp;
            }
            w = w + w * wp;
        }
        mmax <<= 1;
    }
    return;

}

double max_fft_frequency(double _Complex* signal, size_t N)
{
    size_t M = N;
    //double _Complex* fft_spectrum = malloc(N*sizeof(double _Complex));
    double _Complex* fft_spectrum = init_fft(signal, N, &M);
    fft2(fft_spectrum, M);

    size_t imax = 0;
    double max = 0;

    for(size_t i = 0; i < M/2; i++)
    {
        double amp = creal(fft_spectrum[i])*creal(fft_spectrum[i]) + cimag(fft_spectrum[i])*cimag(fft_spectrum[i]); 
        if( amp > max )
        {
            max = amp;
            imax = i;
        } 
    }

    free(fft_spectrum);
    return (1.0*imax)/M;
}

double max_fft_frequency_all(double _Complex* signal, size_t N)
{
    size_t M = N;
    //double _Complex* fft_spectrum = malloc(N*sizeof(double _Complex));
    double _Complex* fft_spectrum = init_fft(signal, N, &M);
    fft2(fft_spectrum, M);

    int imax = 0;
    double max = 0;

    for(size_t i = 0; i < M; i++)
    {
        double amp = creal(fft_spectrum[i])*creal(fft_spectrum[i]) + cimag(fft_spectrum[i])*cimag(fft_spectrum[i]); 
        if( amp > max )
        {
            max = amp;
            imax = (int)i;
        } 
    }
    if( imax > M/2 )
        imax -= M;

    free(fft_spectrum);
    return (1.0*imax)/M;
}

double _Complex* init_fft(double _Complex* signal, size_t N, size_t *M)
{
    // round to next power of two
    *M = *M - 1;
    *M |= *M >> 1;
    *M |= *M >> 2;
    *M |= *M >> 4;
    *M |= *M >> 8;
    *M |= *M >> 16;
    *M = *M + 1;
    //===========================

    size_t bits = bitSize(*M);

    double _Complex* fft_spectrum = (double _Complex*)malloc(*M*sizeof(double _Complex));
    for(size_t i = N; i--;)
        fft_spectrum[reverseBits(i,bits)] = signal[i];
    for(size_t i = N; i < *M; i++)
        fft_spectrum[reverseBits(i,bits)] = 0;

    return fft_spectrum;
}
