#include "fft.h"

void separate(double _Complex* A, size_t N)
{
    size_t N2 = N/2;
    double _Complex* B = (double _Complex*)malloc(N2*sizeof(double _Complex));
    for( size_t i = N2; i--; )
        B[i] = A[2*i+1];
    for( size_t i = N2; i--; )
        A[i] = A[2*i];
    for( size_t i = N2; i--; )
        A[i+N2] = B[i];
    free(B);
    return;
}

void fft2(double _Complex* X, size_t N)
{
    if( N>=2 )
    {
        double _Complex omega = -I*((2*pi)/N);
        size_t N2 = N/2;
        separate(X,N);
        fft2(X, N/2);
        fft2(X+N2, N/2);
        for( size_t i = N2; i--; )
        {
            double _Complex e = X[i];
            double _Complex o = X[i+N2];
            double _Complex w = cexp(omega*i);
            X[i] = e + w * o;
            X[i+N2] = e - w * o;
        }
    }
}

double max_fft_frequency(double _Complex* signal, size_t N)
{
    size_t M = N;
    //double _Complex* fft_spectrum = malloc(N*sizeof(double _Complex));
    double _Complex* fft_spectrum = init_fft(signal, N, &M);
    fft2(fft_spectrum, N);

    size_t imax = 0;
    double max = 0;

    for(size_t i = N/2; i--;)
    {
        if( cabs(fft_spectrum[i]) > max )
        {
            max = cabs(fft_spectrum[i]);
            imax = i;
        } 
    }

    free(fft_spectrum);
    //return (1.0*imax)/N;
    return (1.0*imax)/N;
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

    double _Complex* fft_spectrum = (double _Complex*)malloc(*M*sizeof(double _Complex));
    for(size_t i = N; i--;)
        fft_spectrum[i] = signal[i];
    for(size_t i = N; i < *M; i++)
        fft_spectrum[i] = 0;

    return fft_spectrum;
}

#ifndef COMPILE_WITHOUT_FFTW

double max_fftw_frequency_all(double _Complex* signal, size_t N)
{
    fftw_complex *out = (fftw_complex*) fftw_malloc(N*sizeof(fftw_complex));
    fftw_plan p = fftw_plan_dft_1d(N, signal, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    int imax = 0;
    double max = 0;

    for(size_t i = N; i--;)
    {
        if( cabs(out[i]) > max )
        {
            max = cabs(out[i]);
            imax = (int)i;
        } 
    }
    if( imax > N/2)
        imax = imax - N;

    fftw_destroy_plan(p);
    fftw_free(out);
    return (1.0*imax)/N;
}

double max_fftw_frequency(double _Complex* signal, size_t N)
{
    fftw_complex *out = (fftw_complex*) fftw_malloc(N*sizeof(fftw_complex));
    fftw_plan p = fftw_plan_dft_1d(N, signal, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    size_t imax = 0;
    double max = 0;

    for(size_t i = N/2; i--;)
    {
        if( cabs(out[i]) > max )
        {
            max = cabs(out[i]);
            imax = i;
        } 
    }

    fftw_destroy_plan(p);
    fftw_free(out);
    return (1.0*imax)/N;
}

#endif
