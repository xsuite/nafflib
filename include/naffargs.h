#ifndef __NAFFLIB_NAFFARGS_H__
#define __NAFFLIB_NAFFARGS_H__

struct merit_args
{
    size_t N;
    double _Complex* window;
    double _Complex* signal;
};

typedef struct merit_args merit_args;

#endif
