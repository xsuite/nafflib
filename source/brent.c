#include "brent.h"

double brent_minimize(double (*f)(double,const merit_args*), double min, double max, const merit_args* S)
{

    const int max_iter = 10000;
    const double golden = 0.3819660;
    const double tolerance = 1.490116e-8; //ldexp(1.-25)

    double x, w, v, u; 
    double fu, fv, fw, fx;
    double mid;
    double delta1, delta2;  
    double tol0 = tolerance*0.25;
    double tol1, tol2;  // minimal relative movement in x


    x = w = v = max;
    fw = fv = fx = (*f)(x,S);
    delta2 = delta1 = 0;


    for(int i = max_iter; i--;)
    {
        mid = 0.5 * (min + max);
        tol1 = tolerance * fabs(x) + tol0;
        tol2 = 2. * tol1;
        if( fabs(x - mid) <= (tol2 - 0.5 * (max - min)) )
            return x;

        if( fabs(delta2) > tol1 )
        {
            double r = (x - w) * (fx - fv);
            double q = (x - v) * (fx - fw);
            double p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if( q > 0 )
                p = -p;
            q = fabs(q);
            double delta0 = delta2;
            delta2 = delta1;
            if( fabs(p) >= fabs((0.5 * q) * delta0) || p <= q * (min - x) || p >= q * (max - x) )
            {
                delta2 = (x >= mid) ? min - x : max - x;
                delta1 = golden * delta2;
            }
            else
            {   
                delta1 = p / q;
                u = x + delta1;
                if( (u - min) < tol2 || (max - u) < tol2)
                    delta1 = (mid - x) < 0 ? -fabs(tol1) : fabs(tol1);
            }
        }
        else
        {
            delta2 = (x >= mid) ? min - x : max - x;
            delta1 = golden * delta2;
        }
        u = (fabs(delta1) >= tol1) ? x + delta1 : x + ( (delta1 > 0) ? fabs(tol1) : -fabs(tol1) );
        fu = (*f)(u,S);
        if(fu <= fx)
        {
            if( u >= x )
                min = x;
            else
                max = x;
            v = w;
            w = x;
            x = u;
            fv = fw;
            fw = fx;
            fx = fu;
        }
        else
        {
            if( u < x )
                min = u;
            else
                max = u;
            if( fu <= fw || w == x )
            {
                v = w;
                w = u;
                fv = fw;
                fw = fu;
            }
            else if( fu <= fv || v == x || v == w)
            {
                v = u;
                fv = fu;
            }
        }
    }
    printf("WARNING: nafflib Brent minimization reached maximum number of iterations: %d.\n",max_iter);
    return x;
}

