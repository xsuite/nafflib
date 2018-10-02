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

void py_f1(double _Complex* signal, size_t N, double order, double fft_estimate, double* tune)
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
    *tune = naff_estimate;
    return;
}

void pyget_f1(double _Complex* signal, size_t N, double order, double *tune)
{
    merit_args *margs = (merit_args*)malloc(sizeof(merit_args));
    margs->N = N;
    margs->window = (double _Complex*)malloc(N*sizeof(double _Complex));
    margs->signal = signal;
    hann_harm_window(margs->window, N, order);
    strip_DC(signal, N);

    double (*merit_function)(double,const merit_args*) = minus_magnitude_fourier_integral;

    double step = 1./N;
    double fft_estimate = max_fftw_frequency(margs->signal, N);
    double naff_estimate = brent_minimize( merit_function, fft_estimate-step, fft_estimate+step, margs); 
    
    free(margs->window);
    free(margs);
    //return naff_estimate;
    *tune = naff_estimate;
    return;
}

void hi()
{
    size_t N = 1001;
    //size_t N = 10001;
    merit_args *margs = (merit_args*)malloc(sizeof(merit_args));
    margs->window = (double _Complex*)malloc(N*sizeof(double _Complex));
    margs->signal = (double _Complex*)malloc(N*sizeof(double _Complex));
    margs->N = N;
    hann_harm_window(margs->window, N, 3);
    for(size_t i =0; i < N; i++) margs->signal[i] = cos(2*pi*0.12345678*i);

    double (*merit_function)(double,const merit_args*) = minus_magnitude_fourier_integral;

    double step = 1./N;
    double first_estimate = max_fft_frequency(margs->signal, N);
    double naff_estimate = brent_minimize( merit_function, first_estimate-step, first_estimate+step, margs);  
    while(fabs(fabs(naff_estimate - first_estimate) - step) < 1.e-16)
    {
        double sgn = naff_estimate < first_estimate ? -1 : +1;
        first_estimate += sgn*step; 
        naff_estimate = brent_minimize( merit_function, first_estimate-step, first_estimate+step, margs);  
        printf("edge encounter %.16lf  %.5lf\n",fabs(naff_estimate-first_estimate)-step,sgn);
    }
    
    printf("%lf \n", first_estimate);
    printf("%.10lf \n", naff_estimate);
//    printf(" x[5] %.10lf \n", creal(margs->signal[5]));

    free(margs->signal);
    free(margs->window);
    free(margs);
    return;
}
