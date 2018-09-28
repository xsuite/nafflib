#include "frequency.h"

double minus_magnitude_fourier_integral(double frequency, const merit_args* S)
{
    return -cabs(inner_product(S->signal, 1., frequency, S->window, S->N));
}

double get_f1(double _Complex* signal, size_t N, double order, double fft_estimate)
{
    merit_args *margs = (merit_args*)malloc(sizeof(merit_args));
    margs->N = N;
    margs->window = (double _Complex*)malloc(N*sizeof(double _Complex));
    margs->signal = signal;
    hann_harm_window(margs->window, N, order);
    strip_DC(signal, N);

    double (*merit_function)(double,const merit_args*) = minus_magnitude_fourier_integral;

    double step = 1./N;
    double naff_estimate = brent_minimize( merit_function, fft_estimate-step, fft_estimate+step, margs); 
    
//    printf("%lf \n", fft_estimate);
//    printf("%.10lf \n", naff_estimate);
    free(margs->window);
    free(margs);
    return naff_estimate;
}

void pyget_f1(double _Complex* signal, size_t N, double order, double *tune)
{
    merit_args *margs = (merit_args*)malloc(sizeof(merit_args));
    margs->N = N;
    margs->window = (double _Complex*)malloc(N*sizeof(double _Complex));
    margs->signal = signal;
    hann_harm_window(margs->window, N, order);
    //strip_DC(signal, N);

    double (*merit_function)(double,const merit_args*) = minus_magnitude_fourier_integral;

    double step = 1./N;
    double fft_estimate = max_fft_frequency(margs->signal, N);
    double naff_estimate = brent_minimize( merit_function, fft_estimate-step, fft_estimate+step, margs); 
    
    printf("%.10lf\n",naff_estimate);
//    printf("%lf \n", fft_estimate);
//    printf("%.10lf \n", naff_estimate);
    free(margs->window);
    free(margs);
    //return naff_estimate;
    *tune = naff_estimate;
    return;
}

void hi()
{
    size_t N = 2049;
    merit_args *margs = (merit_args*)malloc(sizeof(merit_args));
    margs->window = (double _Complex*)malloc(N*sizeof(double _Complex));
    margs->signal = (double _Complex*)malloc(N*sizeof(double _Complex));
    margs->N = N;
    hann_harm_window(margs->window, N, 3);
    for(size_t i =0; i < N; i++) margs->signal[i] = cos(2*pi*0.12345678*i)+0.1*cos(2*pi*0.2341*i);

    double (*merit_function)(double,const merit_args*) = minus_magnitude_fourier_integral;

    double step = 1./N;
    double fft_estimate = max_fft_frequency(margs->signal, N);
    double naff_estimate = brent_minimize( merit_function, fft_estimate-step, fft_estimate+step, margs); 
    
    printf("%lf \n", fft_estimate);
    printf("%.10lf \n", naff_estimate);

    free(margs->signal);
    free(margs->window);
    free(margs);
    return;
}
