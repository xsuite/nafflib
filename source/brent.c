#include "brent.h"

double brent_minimize(double (*f)(double,const merit_args*), double min, double max, const merit_args* S)
{

    const int max_iter = 10000;
    const double golden = 0.3819660;
    const double tolerance = 1.490116e-8; //ldexp(1.-25)

    double x;  // minima so far
    double w;  // second best point
    double v;  // previous value of w
    double u;  // most recent evaluation point
    double delta;  // The distance moved in the last step
    double delta2; // The distance moved in the step before last
    double fu, fv, fw, fx;  // function evaluations at u, v, w, x
    double mid; // midpoint of min and max
    double fract1, fract2;  // minimal relative movement in x

    x = w = v = max;
    fw = fv = fx = (*f)(x,S);
    delta2 = delta = 0;


    for(int i = max_iter; i--;)
    {
        // get midpoint
        mid = (min + max) / 2;
        // work out if we're done already:
        fract1 = tolerance * fabs(x) + tolerance / 4;
        fract2 = 2 * fract1;
        if(fabs(x - mid) <= (fract2 - (max - min) / 2))
            return x;

        if(fabs(delta2) > fract1)
        {
            // try and construct a parabolic fit:
            double r = (x - w) * (fx - fv);
            double q = (x - v) * (fx - fw);
            double p = (x - v) * q - (x - w) * r;
            q = 2 * (q - r);
            if(q > 0)
                p = -p;
            q = fabs(q);
            double td = delta2;
            delta2 = delta;
            // determine whether a parabolic step is acceptible or not:
            if((fabs(p) >= fabs(q * (td / 2))) || (p <= q * (min - x)) || (p >= q * (max - x)))
            {
                // nope, try golden section instead
                delta2 = (x >= mid) ? min - x : max - x;
                delta = golden * delta2;
            }
            else
            {
                // whew, parabolic fit:
                delta = p / q;
                u = x + delta;
                if(((u - min) < fract2) || ((max- u) < fract2))
                    delta = (mid - x) < 0 ? -fabs(fract1) : fabs(fract1);
            }
        }
        else
        {
            // golden section:
            delta2 = (x >= mid) ? min - x : max - x;
            delta = golden * delta2;
        }
        // update current position:
        u = (fabs(delta) >= fract1) ? x + delta : (delta > 0 ? x + fabs(fract1) : x - fabs(fract1));
        fu = (*f)(u,S);
        if(fu <= fx)
        {
            // good new point is an improvement!
            // update brackets:
            if(u >= x)
                min = x;
            else
                max = x;
            // update control points:
            v = w;
            w = x;
            x = u;
            fv = fw;
            fw = fx;
            fx = fu;
        }
        else
        {
            // Oh dear, point u is worse than what we have already,
            // even so it *must* be better than one of our endpoints:
            if(u < x)
                min = u;
            else
                max = u;
            if((fu <= fw) || (w == x))
            {
                // however it is at least second best:
                v = w;
                w = u;
                fv = fw;
                fw = fu;
            }
            else if((fu <= fv) || (v == x) || (v == w))
            {
                // third best:
                v = u;
                fv = fu;
            }
        }

    }
    printf("WARNING: nafflib Brent minimization reached maximum number of iterations: %d.\n",max_iter);
    return x;
}
