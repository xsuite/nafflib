#include "signal_processing.h"

double _Complex inner_product(const double _Complex* signal, double amplitude, double frequency,const double _Complex* window, size_t N)
{
    double omega = (2*pi)*frequency;
    double _Complex result = 0.;
    double theta = omega*N;
    for( size_t i = N; i--; )
    {
        theta -= omega;
        result += (signal[i]*(cos(theta)-I*sin(theta)))*window[i];
    }
    return (amplitude*result)/N;
}

void strip_DC(double _Complex* signal, size_t N)
{
    double _Complex mean = 0.;
    for(int i = N; i--;)
        mean += signal[i];
    mean /= N;
    for(int i = N; i--;)
        signal[i] -= mean;
    return;
}

double _Complex frequency_project(double f0, double _Complex *amps, double *freqs, size_t n_frequencies, double _Complex *window, size_t N)
{
    double _Complex num = 0.;
    double _Complex den = 0.;

    double omega0 = (2*pi)*f0;
    double omega[n_frequencies];
    for(size_t j = n_frequencies; j--;)
        omega[j] = (2*pi)*freqs[j];

    for(size_t i = N; i--;)
    {
        double _Complex Ai = 0;
        for(size_t j = n_frequencies; j--;)
         //   Ai += amps[j]*cexp(I*(omega[j]*i));
            Ai += amps[j]*(cos(omega[j]*i)+I*sin(omega[j]*i));
        double _Complex cAi = creal(Ai) - I*cimag(Ai);

        //num += (cexp(I*omega0*i) * cAi) * window[i]; 
        num += ((cos(omega0*i)+I*sin(omega0*i)) * cAi) * window[i]; 
        den += (Ai * cAi) * window[i]; 
    }

    return (num/den);
}

void subtract_frequency(double _Complex* signal, double _Complex** normal_amps, double* frequencies, double _Complex* amplitudes, size_t i, double _Complex* window, size_t N)
{
    for(size_t k = 0; k < i + 1; k++)
        normal_amps[i][k] = 0.;
    for(size_t j = 0; j < i; j++)
    {
        double Aij = frequency_project(frequencies[i], normal_amps[j], frequencies, j+1, window, N);
        for(size_t k = 0; k < j + 1; k++)
            normal_amps[i][k] -= Aij*normal_amps[j][k];
    }
    normal_amps[i][i]=1;
    double _Complex frequency_amplitude = signal_project(signal, normal_amps[i], frequencies, i+1, window, N);
    for(size_t k = 0; k < i+1; k++)
        normal_amps[i][k] *= frequency_amplitude;

    remove_component(signal, normal_amps[i], frequencies, i+1, N); 
    amplitudes[i] = frequency_amplitude;
    return;

}

void remove_component( double _Complex* signal, double _Complex* amps, double* frequencies, size_t n_freqs, size_t N)
{
    double omega[n_freqs];
    for(size_t j = n_freqs; j--;)
        omega[j] = (2*pi)*frequencies[j];

    for(size_t i = N; i--;)
    {
        for(size_t j = n_freqs; j--;)
            signal[i] -= amps[j]*(cos(omega[j]*i)+I*sin(omega[j]*i)); 
    }
    return;
}

double _Complex signal_project(double _Complex* signal, double _Complex* amps, double* freqs, size_t n_frequencies, double _Complex* window, size_t N)
{

    double _Complex num = 0.;
    double _Complex den = 0.;

    double omega[n_frequencies];
    for(size_t j = n_frequencies; j--;)
        omega[j] = (2*pi)*freqs[j];

    for(size_t i = N; i--;)
    {
        double _Complex Ai = 0;
        for(size_t j = n_frequencies; j--;)
            //Ai += amps[j]*cexp(I*(omega[j]*i));
            Ai += amps[j]*(cos(omega[j]*i)+I*sin(omega[j]*i));
        double _Complex cAi = creal(Ai) - I*cimag(Ai);

        num += (signal[i] * cAi) * window[i]; 
        den += (Ai * cAi) * window[i]; 
    }

    return (num/den);

}
