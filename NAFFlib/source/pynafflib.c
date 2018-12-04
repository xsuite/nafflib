#include "pynafflib.h"

#ifndef COMPILE_WITHOUT_PYTHON

static PyObject* get_tune(PyObject* self, PyObject* args)
{
    double order = 2; // default value
    int interpolate_integral = 0;  // default value

    PyObject *arg1 = NULL;
    PyObject *in_array = NULL;
    if(!PyArg_ParseTuple(args, "O!|di", &PyArray_Type, &arg1, &order, &interpolate_integral))
        return NULL;

    in_array = PyArray_FROM_OTF(arg1, NPY_COMPLEX128, NPY_ARRAY_IN_ARRAY); 
    if(in_array == NULL)
        return NULL;
    
    double _Complex* signal = (double _Complex *)PyArray_DATA((PyArrayObject *)in_array);
    int nd = (int)PyArray_NDIM((PyArrayObject *)in_array);
    if(nd != 1)
    {
        printf("***ERROR***: input array is not 1D.\n");
        Py_INCREF(Py_None);
        return Py_None;
    }
    size_t *dims = (size_t *)PyArray_DIMS((PyArrayObject *)in_array);

    double tune = get_q(signal, dims[0], order, interpolate_integral);

    Py_DECREF(in_array);
    return Py_BuildValue("d", tune);
}

static PyObject* get_tunes(PyObject* self, PyObject* args)
{
    double order = 2; // default value
    int interpolate_integral = 0;  // default value
    int nfreqs = 1;  // default value

    PyObject *arg1 = NULL;
    PyObject *in_array = NULL;
    if(!PyArg_ParseTuple(args, "O!|idi", &PyArray_Type, &arg1, &nfreqs, &order, &interpolate_integral))
        return NULL;

    in_array = PyArray_FROM_OTF(arg1, NPY_COMPLEX128, NPY_ARRAY_IN_ARRAY); 
    if(in_array == NULL)
        return NULL;
    
    int nd = (int)PyArray_NDIM((PyArrayObject *)in_array);
    if(nd != 1)
    {
        printf("***ERROR***: input array is not 1D.\n");
        Py_INCREF(Py_None);
        return Py_None;
    }
    size_t *dims = (size_t *)PyArray_DIMS((PyArrayObject *)in_array);
    if(dims[0] == 0)
    {
        printf("***ERROR***: input array is empty.\n");
        Py_INCREF(Py_None);
        return Py_None;
    }

    double _Complex* in_signal = (double _Complex *)PyArray_DATA((PyArrayObject *)in_array);
    double _Complex* signal = (double _Complex *)malloc(dims[0]*sizeof(double _Complex));
    for(int i = dims[0]; i--;)
        signal[i] = in_signal[i];

    double *frequencies = (double *)malloc(nfreqs*sizeof(double));
    double _Complex *positive_amplitudes = (double _Complex*)malloc(nfreqs*sizeof(double _Complex));
    double _Complex *negative_amplitudes = (double _Complex*)malloc(nfreqs*sizeof(double _Complex));

    get_f_neg(signal, dims[0], order, frequencies, positive_amplitudes, negative_amplitudes, nfreqs);
    free(signal);

    npy_intp new_dims[1];
    new_dims[0] = nfreqs;

    PyObject* freqs = PyArray_SimpleNewFromData(1, new_dims, NPY_FLOAT64, (void*) frequencies);
    PyObject* pos_amps = PyArray_SimpleNewFromData(1, new_dims, NPY_COMPLEX128, (void*) positive_amplitudes);
    PyObject* neg_amps = PyArray_SimpleNewFromData(1, new_dims, NPY_COMPLEX128, (void*) negative_amplitudes);

    PyArray_ENABLEFLAGS((PyArrayObject *) freqs, NPY_ARRAY_OWNDATA);
    PyArray_ENABLEFLAGS((PyArrayObject *) pos_amps, NPY_ARRAY_OWNDATA);
    PyArray_ENABLEFLAGS((PyArrayObject *) neg_amps, NPY_ARRAY_OWNDATA);

    Py_DECREF(in_array);
    return Py_BuildValue("NNN", freqs, pos_amps, neg_amps);
}

static PyObject* get_tunes_all(PyObject* self, PyObject* args)
{
    double order = 2; // default value
    int interpolate_integral = 0;  // default value
    int nfreqs = 1;  // default value

    PyObject *arg1 = NULL;
    PyObject *in_array = NULL;
    if(!PyArg_ParseTuple(args, "O!|idi", &PyArray_Type, &arg1, &nfreqs, &order, &interpolate_integral))
        return NULL;

    in_array = PyArray_FROM_OTF(arg1, NPY_COMPLEX128, NPY_ARRAY_IN_ARRAY); 
    if(in_array == NULL)
        return NULL;
    
    int nd = (int)PyArray_NDIM((PyArrayObject *)in_array);
    if(nd != 1)
    {
        printf("***ERROR***: input array is not 1D.\n");
        Py_INCREF(Py_None);
        return Py_None;
    }
    size_t *dims = (size_t *)PyArray_DIMS((PyArrayObject *)in_array);
    if(dims[0] == 0)
    {
        printf("***ERROR***: input array is empty.\n");
        Py_INCREF(Py_None);
        return Py_None;
    }

    double _Complex* in_signal = (double _Complex *)PyArray_DATA((PyArrayObject *)in_array);
    double _Complex* signal = (double _Complex *)malloc(dims[0]*sizeof(double _Complex));
    for(int i = dims[0]; i--;)
        signal[i] = in_signal[i];

    double *frequencies = (double *)malloc(nfreqs*sizeof(double));
    double _Complex *amplitudes = (double _Complex*)malloc(nfreqs*sizeof(double _Complex));

    get_f_all(signal, dims[0], order, frequencies, amplitudes, nfreqs);
    free(signal);

    npy_intp new_dims[1];
    new_dims[0] = nfreqs;

    PyObject* freqs = PyArray_SimpleNewFromData(1, new_dims, NPY_FLOAT64, (void*) frequencies);
    PyObject* amps = PyArray_SimpleNewFromData(1, new_dims, NPY_COMPLEX128, (void*) amplitudes);

    PyArray_ENABLEFLAGS((PyArrayObject *) freqs, NPY_ARRAY_OWNDATA);
    PyArray_ENABLEFLAGS((PyArrayObject *) amps, NPY_ARRAY_OWNDATA);

    Py_DECREF(in_array);
    return Py_BuildValue("NN", freqs, amps);
}

static PyMethodDef NAFFlibMethods[] =
{
    {"get_tune", get_tune, METH_VARARGS, "returns single (most dominant) frequency"},
    {"get_tunes", get_tunes, METH_VARARGS, "returns N dominant frequencies with their amplitudes as well as the amplitudes of the opposite frequencies. Recommended with real-valued inputs."},
    {"get_tunes_all", get_tunes_all, METH_VARARGS, "returns N dominant frequencies with their amplitudes as well, treating negative and positive frequencies separately. Should be used with a complex-valued input."},
    {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef cModPyNAFFlib_c = 
{
    PyModuleDef_HEAD_INIT,
    "NAFFlib_c","The NAFF algorithm written in a python-wrapped C library.",
    -1,
    NAFFlibMethods
};

PyMODINIT_FUNC PyInit_NAFFlib_c(void)
{
    import_array();
    return PyModule_Create(&cModPyNAFFlib_c);
}

#else

PyMODINIT_FUNC initNAFFlib2_c(void)
{
    (void) Py_InitModule("NAFFlib2_c", NAFFlibMethods);
    import_array();
}


#endif

#endif
