// interface between Python and the C++ implementation
// https://docs.microsoft.com/de-de/visualstudio/python/working-with-c-cpp-python-in-visual-studio?view=vs-2019
// https://docs.scipy.org/doc/numpy-1.13.0/user/c-info.how-to-extend.html

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#ifdef _DEBUG
	#undef _DEBUG
	#include <Python.h>
	#define _DEBUG
#else
	#include <Python.h>
#endif

#include <exception>
#include <cstdio>
#include "numpy/noprefix.h"
#include "factoredmarcher.h"


static PyObject* fast_marching_(PyObject* self, PyObject* args, const bool factored)
{
	// define placeholders
	PyObject* pc, * px_s, * pdx;
	int order;

	// parse arguments
	if (!PyArg_ParseTuple(args, "OOOi", &pc, &px_s, &pdx, &order))
		return NULL;

	// only orders 1 and 2 are implemented
	if (order < 1 || order > 2)
	{
		PyErr_SetString(PyExc_ValueError, "currently only orders 1 and 2 are supported.");
		return NULL;
	}

	// convert python-numpy-array to c-numpy-array
	PyArrayObject* c = (PyArrayObject*)PyArray_FROMANY(pc, NPY_DOUBLE, 1, 0, NPY_ARRAY_IN_ARRAY);
	if (!c)
	{
		PyErr_SetString(PyExc_ValueError, "c must be an array of doubles.");
		return NULL;
	}

	int ndim = PyArray_NDIM(c);

	// check if px0 is a python-numpy-array and if it has the correct size
	PyArrayObject* x_s_ = (PyArrayObject*)PyArray_FROMANY(px_s, NPY_LONG, 1, 1, NPY_ARRAY_IN_ARRAY);
	if (!x_s_)
	{
		PyErr_SetString(PyExc_ValueError, "x_s must be a 1D array of ints.");
		Py_DECREF(c);
		return NULL;
	}
	else if (PyArray_SIZE(x_s_) != ndim)
	{
		char msg[100];
		std::sprintf(msg, "size of x_s and number of dimensions of c do not match: %d != %d.", PyArray_SIZE(x_s_), ndim);
		PyErr_SetString(PyExc_ValueError, msg);
		Py_DECREF(c);
		return NULL;
	}

	// check if pdx is a python-numpy-array and if it has the correct size
	PyArrayObject* dx_ = (PyArrayObject*)PyArray_FROMANY(pdx, NPY_DOUBLE, 1, 1, NPY_ARRAY_IN_ARRAY);
	if (!dx_)
	{
		PyErr_SetString(PyExc_ValueError, "dx must be a 1D array of doubles.");
		Py_DECREF(c);
		Py_DECREF(x_s_);
		return NULL;
	}
	else if (PyArray_SIZE(dx_) != ndim)
	{
		char msg[100];
		std::sprintf(msg, "size of dx and number of dimensions of c do not match: %d != %d.", PyArray_SIZE(dx_), ndim);
		PyErr_SetString(PyExc_ValueError, msg);
		Py_DECREF(c);
		return NULL;
	}

	double* dx = (double*)PyArray_DATA(dx_);
	long* shape = new long[ndim];
	// index version of x0_
	size_t x_s = 0;

	long* x_s_d = (long*)PyArray_DATA(x_s_);
	size_t tmp = 1;
	for (int i = ndim - 1; i >= 0; i--)
	{
		shape[i] = (long)PyArray_DIM(c, i);

		if (x_s_d[i] < 0 || x_s_d[i] >= shape[i])
		{
			char msg[128];
			std::sprintf(msg, "entries of x_s have to be within the shape of c, but entry %d with value %d is not in [0, %d].", i, x_s_d[i], shape[i]-1);
			PyErr_SetString(PyExc_ValueError, msg);
			Py_DECREF(c);
			Py_DECREF(x_s_);
			Py_DECREF(dx_);
			delete[] shape;
			return NULL;
		}
		x_s += tmp * x_s_d[i];
		tmp *= shape[i];

		if (dx[i] <= 0)
		{
			char msg[128];
			std::sprintf(msg, "entries of dx must be greater than zero, but entry %d with value %f is not.", i, dx[i]);
			PyErr_SetString(PyExc_ValueError, msg);
			Py_DECREF(c);
			Py_DECREF(x_s_);
			Py_DECREF(dx_);
			delete[] shape;
			return NULL;
		}
	}

	// create a new array for the return value
	PyArrayObject* tau = (PyArrayObject*)PyArray_ZEROS(ndim, PyArray_DIMS(c), NPY_DOUBLE, 0);
	if (!tau)
	{
		Py_DECREF(c);
		Py_DECREF(x_s_);
		Py_DECREF(dx_);
		delete[] shape;
		return NULL;
	}

	Marcher* m;
	if (factored)
		m = new FactoredMarcher((double*)PyArray_DATA(c), ndim, shape, dx, order);
	else
		m = new Marcher((double*)PyArray_DATA(c), ndim, shape, dx, order);


	try
	{
		m->solve(x_s, (double*)PyArray_DATA(tau));
	}
	catch (const std::exception& ex)
	{
		// propagate error
		PyErr_SetString(PyExc_RuntimeError, ex.what());
		Py_XDECREF(c);
		Py_XDECREF(x_s_);
		Py_XDECREF(dx_);
		delete m;
		delete[] shape;

		return NULL;
	}

	delete m;
	delete[] shape;

	return (PyObject*)tau;
}

static PyObject* fast_marching_wrapper(PyObject* self, PyObject* args)
{
	return fast_marching_(self, args, false);
}

static PyObject* factored_marching_wrapper(PyObject* self, PyObject* args)
{
	return fast_marching_(self, args, true);
}


// Each entry in this array is a PyMethodDef structure containing 
// 1) the Python name,
// 2) the C - function that implements the function,
// 3) flags indicating whether or not keywords are accepted for this function, and 
// 4) The docstring for the function.
static PyMethodDef fm_methods[] = {
	{"fast_marching", (PyCFunction)fast_marching_wrapper, METH_VARARGS, NULL},
	{"factored_fast_marching", (PyCFunction)factored_marching_wrapper, METH_VARARGS, NULL},

	// Terminate the array with an object containing nulls.
	{NULL, NULL, 0, NULL}
};


static PyModuleDef cfm_module = {
	PyModuleDef_HEAD_INIT,
	"cfm",																		// Module name to use with Python import statements
	"c++ implementation of (factored) fast marching for the eikonal equation",  // Module description
	0,
	fm_methods																	// Structure that defines the methods of the module
};


PyMODINIT_FUNC PyInit_cfm(void)
{
	PyObject* m = PyModule_Create(&cfm_module);
	if (m == NULL)
		return NULL;

	import_array();
	return m;
}