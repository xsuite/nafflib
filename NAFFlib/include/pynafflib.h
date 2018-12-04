#ifndef __NAFFLIB_PYNAFFLIB_H__
#define __NAFFLIB_PYNAFFLIB_H__

#ifndef COMPILE_WITHOUT_PYTHON

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>
#include <complex.h>

#include "frequency.h"

static PyObject* get_tune(PyObject* self, PyObject* args);
static PyObject* get_tunes(PyObject* self, PyObject* args);
static PyObject* get_tunes_all(PyObject* self, PyObject* args);

#endif

#endif
