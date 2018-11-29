#include "windows.h"

double cheb_poly(int n, double x)
{
    double res;

    if(fabs(x) <= 1)
        res = cos(n*acos(x));
    else
        res = cosh(n*acosh(x));

    return res;
}

void cheb_window(double _Complex window[], const size_t N, const double a){
    assert( window != NULL );
    double max = 0;
    double tg = pow(10, a);
    double x0 = cosh((1.0/(N-1))*acosh(tg));
    double M = (N-1)/2;
    
    if( N%2 == 0 )
        M += 0.5;

    for( size_t nn = 0; nn < (N/2 + 1); nn++)
    {
        double n = nn - M;
        double sum = 0;

        for( size_t i = 1; i <= M; i++ )
            sum += cheb_poly(N-1, x0*cos((pi*i)/N))*cos((2.0*n)*((pi*i)/N));

        window[nn] = tg + 2*sum;
        window[N-nn-1] = window[nn];
        
        if( creal(window[nn]) > max )
            max = creal(window[nn]);
    }
    for(size_t nn = 0; nn < N; nn++)
        window[nn] = window[nn]/max;

    return;
}

void hann_harm_window(double _Complex window[], const size_t N, const double n)
{
    double T1 = 0.;
    double T2 = N;
    double TM = (T2-T1)/2.;
    double PIST = pi/TM;

    int factorial_1 = 1;
    int factorial_2 = 1;

    for( size_t j = 1; j <= n; j++ )
        factorial_1 *= j;
    for( size_t j = 1; j <= 2*n; j++ )
        factorial_2 *= j;

    //double cn = pow(2., n)*(((1.*factorial_1)*factorial_1)/(1.*factorial_2));
    double cn = exp(n*log(2.))*(((1.*factorial_1)*factorial_1)/(1.*factorial_2));
    
    for( size_t i = 0; i < N; i++)
        //window[i] = cn*pow(1. + cos( (i-TM)*PIST ) ,n);
        window[i] = cn*exp(n*log(1. + cos( (i-TM)*PIST )));
    
    return;
}


double calculateFm(size_t m, double sp2, double a, size_t nBar)
{
    int n[nBar-1];
    int p[nBar-1];

    double Num[nBar-1];
    double Den[nBar-1];
//  double Fm[nBar-1];

    double prodNum = 1.;
    double prodDen = 1.;

    for( size_t i = 0; i < nBar - 1; i++ )
    {
        n[i] = i + 1;
        p[i] = i+1;
//        Fm[i] = 0;
        Num[i] = 1 -  ((m*m)/sp2)/( a * a + (n[i] - 0.5) * (n[i] - 0.5) );
        prodNum *= Num[i];
        Den[i] = 1. -  (m*m)/((1.*p[i])*p[i]);
        if( (i+1) != m )
            prodDen *= Den[i];
    }
    double prodFm = pow(-1,(m+1))*(prodNum/(2*prodDen));
    return prodFm;
}


void taylorWindow( double _Complex w[], const size_t N, const double NBAR, const double SLL)
{
    double A = acosh( pow(10,SLL) )/pi;
    double SP2 = NBAR*(NBAR/(A*A + (NBAR-0.5)*(NBAR-0.5)));
    double Fm[(int)NBAR], Xi[N];
    if( NBAR < (2 * pow(A,2) + 0.5) ) {
        printf("Warning (taylorWindow): SSL value does nto satisfy NBAR >= 2*A^2+0.5\n");
        printf(" A = %lf\n", A);
        printf(" SLL = %lf\n", SLL);
        printf(" NBAR = %lf\n", NBAR);
        printf(" NBAR should be >= %lf\n", 2.*pow(A,2)+0.5);
    }
    
    for(size_t i = 0; i < N; i++)
    {
        w[i] = 1.;
        Xi[i] = (i - 0.5*N + 0.5)/N;
    }

    for( size_t i = 0; i < NBAR - 1; i++)
        Fm[i] = calculateFm( i+1, SP2, A, NBAR );

    for( size_t i = 0; i < N; i++)
    {
        double sum = 0.;
        for( size_t j = 1; j < NBAR; j++ )
            sum += Fm[j-1] * cos( ((2*pi)*j)*Xi[i] );
        w[i] += 2*sum;
    }
    return;
}



void taylor_window(double _Complex window[], const size_t N, const double SLL)
{
    double A = acosh( pow(10,SLL) )/pi;
    double NBAR = round( (2 * pow(A,2) + 0.5) + 0.5 );

    taylorWindow( window, N, NBAR, SLL);
    return;
}


void no_window( double _Complex window[], const size_t N )
{
    for( size_t i = 0; i < N; i++)
        window[i] = 1.;
    return;
}












