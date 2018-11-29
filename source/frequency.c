#include "frequency.h"

double minus_magnitude_fourier_integral(double frequency, const merit_args* S)
{
    double _Complex amp = inner_product(S->signal, 1., frequency, S->window, S->N);
    return -(creal(amp)*creal(amp)+cimag(amp)*cimag(amp));
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
    
    free(margs->window);
    free(margs);
    return naff_estimate;
}

size_t interpolating_size(size_t N)
{
    while(N%6 != 1)
        N--;
    return N;
}

void use_interpolating_integral_hardy(size_t N, double _Complex* window)
{
    //Modify window (weighting function) to correspond to a 7-point newton cotes-like interpolation when summing
    double h = 1./N;
    double c1 = (28./100.)*h;
    double c2 = (56./100.)*h;
    double c3 = (162./100.)*h;
    double c4 = (0./140.)*h;
    double c5 = (220./100.)*h;

    window[0] *= c1;
    for(size_t i = 0; i < N-1; i+=6)
    {
        window[i+1] *= c3;
        window[i+2] *= c4;
        window[i+3] *= c5;
        window[i+4] *= c4;
        window[i+5] *= c3;
        window[i+6] *= c2;
    }
    window[N-1] *= 0.5;
    return;
}

void use_interpolating_integral(size_t N, double _Complex* window)
{
    //Modify window (weighting function) to correspond to a 7-point newton cotes interpolation when summing
    double h = 1./N;
    double c1 = (41./140.)*h;
    double c2 = (82./140.)*h;
    double c3 = (216./140.)*h;
    double c4 = (27./140.)*h;
    double c5 = (272./140.)*h;

    window[0] *= c1;
    for(size_t i = 0; i < N-1; i+=6)
    {
        window[i+1] *= c3;
        window[i+2] *= c4;
        window[i+3] *= c5;
        window[i+4] *= c4;
        window[i+5] *= c3;
        window[i+6] *= c2;
    }
    window[N-1] *= 0.5;
    return;
}

double get_q(double _Complex* signal, size_t N, double order, int interpolate_integral)
{
    if(interpolate_integral)
        N = interpolating_size(N);
    merit_args *margs = (merit_args*)malloc(sizeof(merit_args));
    margs->N = N;
    margs->window = (double _Complex*)malloc(margs->N*sizeof(double _Complex));
    margs->signal = signal;
    hann_harm_window(margs->window, margs->N, order);
    strip_DC(margs->signal, margs->N);
    if(interpolate_integral)
        use_interpolating_integral(margs->N, margs->window);

    double (*merit_function)(double,const merit_args*) = minus_magnitude_fourier_integral;

    double step = 1./margs->N;
    double fft_estimate = max_fft_frequency(margs->signal, margs->N);
    double naff_estimate = brent_minimize( merit_function, fft_estimate-step, fft_estimate+step, margs); 
    
    free(margs->window);
    free(margs);
    return naff_estimate;
}

void get_f_all(double _Complex* signal, size_t N, double order, double* frequencies, double _Complex* amplitudes, size_t n_freqs)
{
    merit_args *margs = (merit_args*)malloc(sizeof(merit_args));
    margs->N = N;
    margs->window = (double _Complex*)malloc(N*sizeof(double _Complex));
    margs->signal = signal;
    hann_harm_window(margs->window, N, order);
    strip_DC(signal, N);

    double (*merit_function)(double,const merit_args*) = minus_magnitude_fourier_integral;

    double step = 1./N;

    double _Complex* normal_amplitudes[n_freqs];
    for(size_t i = 0; i < n_freqs; i++)
    {
        normal_amplitudes[i] = (double _Complex*)malloc((i+1)*sizeof(double _Complex));
    }
    for(size_t i = 0; i < n_freqs; i++)
    {
        double fft_estimate = max_fft_frequency_all(margs->signal, N);
        double naff_estimate = brent_minimize( merit_function, fft_estimate-step, fft_estimate+step, margs); 
        frequencies[i] = naff_estimate;
        subtract_frequency(margs->signal, normal_amplitudes, frequencies, amplitudes, i, margs->window, N);
    }
    
    free(margs->window);
    free(margs);
    for(size_t i = 0; i < n_freqs; i++)
    {
       free(normal_amplitudes[i]); 
    }
    return;
}

void get_f_neg(double _Complex* signal, size_t N, double order, double* frequencies, double _Complex* amplitudes, double _Complex* negative_amplitudes, size_t n_freqs)
{
    merit_args *margs = (merit_args*)malloc(sizeof(merit_args));
    margs->N = N;
    margs->window = (double _Complex*)malloc(N*sizeof(double _Complex));
    margs->signal = signal;
    hann_harm_window(margs->window, N, order);
    strip_DC(signal, N);

    double (*merit_function)(double,const merit_args*) = minus_magnitude_fourier_integral;

    double step = 1./N;

    double _Complex* normal_amplitudes[n_freqs];
    double _Complex* normal_amplitudes_neg[n_freqs];
    double frequencies_neg[n_freqs];
    for(size_t i = 0; i < n_freqs; i++)
    {
        normal_amplitudes[i] = (double _Complex*)malloc((i+1)*sizeof(double _Complex));
        normal_amplitudes_neg[i] = (double _Complex*)malloc((i+1)*sizeof(double _Complex));
    }
    //double frequencies[n_freqs];
    //double amplitudes[n_freqs];
    for(size_t i = 0; i < n_freqs; i++)
    {
        double fft_estimate = max_fft_frequency(margs->signal, N);
        //double naff_estimate = brent_minimize_integral( fft_estimate-step, fft_estimate+step, margs); 
        //double neg_naff_estimate = brent_minimize_integral( -naff_estimate-step, -naff_estimate+step, margs); 
        double naff_estimate = brent_minimize( merit_function, fft_estimate-step, fft_estimate+step, margs); 
        double neg_naff_estimate = brent_minimize( merit_function, -naff_estimate-step, -naff_estimate+step, margs); 
        frequencies[i] = naff_estimate;
        frequencies_neg[i] = neg_naff_estimate;
        subtract_frequency(margs->signal, normal_amplitudes, frequencies, amplitudes, i, margs->window, N);
        subtract_frequency(margs->signal, normal_amplitudes_neg, frequencies_neg, negative_amplitudes, i, margs->window, N);

    }
    
//    printf("%lf \n", fft_estimate);
//    printf("%.10lf \n", naff_estimate);
    free(margs->window);
    free(margs);
    for(size_t i = 0; i < n_freqs; i++)
    {
       free(normal_amplitudes[i]); 
       free(normal_amplitudes_neg[i]); 
    }
    return;
    //return naff_estimate;
}

