#include "toy1.h"

void hello(int* N)
{
    for(int i = *N; i--;)
        printf("helllooo\n");
    *N = *N + 2;
    return;
}

void hello2(toy_data *data, const int N)
{
    for(int i = N; i--;)
        printf("%d/%d yo %d\n",i,N,data->a);

    return;

}

void hi3(double array[], const int as)
{
    double sum = 0;
    for(int i = as; i--;)
        sum += array[i];
    printf("Sum is %lf \n", sum);
    return;
}
