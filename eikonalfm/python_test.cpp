// memory leak detection
//#include "debugtests.h"

#include <iostream>
#include <iomanip>
//#include <errcode.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#ifdef _DEBUG
	#undef _DEBUG
	#include <Python.h>
	#define _DEBUG
#else
	#include <Python.h>
#endif

//#include "numpy/noprefix.h"
#include "numpy/arrayobject.h"
#include "marcher.hpp"

using namespace std;

static PyObject* fast_marching(PyObject* self, PyObject* args)
{
	// define placeholders
	PyObject* pc, * px0, * pdx;
	int order;

	// parse arguments
	if (!PyArg_ParseTuple(args, "OOOi", &pc, &px0, &pdx, &order))
		return NULL;

	// convert python-numpy-array to c-numpy-array
	PyArrayObject* c = (PyArrayObject*)PyArray_FROMANY(pc, NPY_DOUBLE, 1, 0, NPY_ARRAY_IN_ARRAY);
	if (!c)
	{
		PyErr_SetString(PyExc_ValueError, "c must be an array of doubles");
		return NULL;
	}

	int ndim = PyArray_NDIM(c);

	// check if px0 is a python-numpy-array and if it has the correct size
	PyArrayObject* x0_ = (PyArrayObject*)PyArray_FROMANY(px0, NPY_LONG, 1, 1, NPY_ARRAY_IN_ARRAY);
	if (!x0_ || PyArray_SIZE(x0_) != ndim)
	{
		PyErr_SetString(PyExc_ValueError, "x0 must be an array of ints with size len(c.shape)");
		Py_DECREF(c);
		return NULL;
	}

	// check if pdx is a python-numpy-array and if it has the correct size
	PyArrayObject* dx_ = (PyArrayObject*)PyArray_FROMANY(pdx, NPY_DOUBLE, 1, 1, NPY_ARRAY_IN_ARRAY);
	if (!dx_ || PyArray_SIZE(dx_) != ndim)
	{
		PyErr_SetString(PyExc_ValueError, "dx must be an array of doubles with size len(c.shape)");
		Py_DECREF(c);
		Py_DECREF(x0_);
		return NULL;
	}

	double* dx = (double*)PyArray_DATA(dx_);
	unsigned long* shape = new unsigned long[ndim];
	// index version of x0_
	size_t x0 = 0;

	long* x0_d = (long*)PyArray_DATA(x0_);
	size_t tmp = 1;
	for (int i = ndim - 1; i >= 0; i--)
	{
		shape[i] = (long)PyArray_DIM(c, i);

		if (x0_d[i] < 0 || x0_d[i] >= shape[i])
		{
			PyErr_SetString(PyExc_IndexError, "values in x0 have to be within the shape of c");
			Py_DECREF(c);
			Py_DECREF(x0_);
			Py_DECREF(dx_);
			delete[] shape;
			return NULL;
		}
		x0 += tmp * x0_d[i];
		tmp *= shape[i];

		if (dx[i] <= 0)
		{
			PyErr_SetString(PyExc_IndexError, "dx must be greater than zero");
			Py_DECREF(c);
			Py_DECREF(x0_);
			Py_DECREF(dx_);
			delete[] shape;
			return NULL;
		}
	}

	delete[] x0_d;

	// create a new array for the return value
	PyArrayObject* tau = (PyArrayObject*)PyArray_ZEROS(ndim, PyArray_DIMS(c), NPY_DOUBLE, 0);
	if (!tau)
	{
		Py_DECREF(c);
		Py_DECREF(x0_);
		Py_DECREF(dx_);
		delete[] shape;
		return NULL;
	}

	auto info = MarcherInfo{ndim, shape};
	Marcher* m = new Marcher((double*)PyArray_DATA(c), info, dx, order);

	m->solve(x0, (double*)PyArray_DATA(tau));

	delete m;
	delete[] shape;

	return (PyObject*)tau;
}


int main()
{
	//Py_SetPath(Py_DecodeLocale("E:\\Anaconda\\Lib;E:\\Anaconda\\Lib\\site-packages;E:\\Anaconda\\DLLs", NULL));
    //Py_SetPath(Py_DecodeLocale("~\\anaconda3\\lib;~\\anaconda3\\lib\\python3.7\\site-packages;", NULL));
    Py_SetPath(Py_DecodeLocale("D:\\Funktionsprogramme\\Anaconda\\lib;D:\\Funktionsprogramme\\Anaconda\\lib\\site-packages;D:\\Funktionsprogramme\\Anaconda\\DLLs", NULL));

	Py_Initialize();
	import_array();


	double c_[] = {
		0.50, 0.85, 1.20, 1.55, 1.90, 2.25, 2.60, 2.95, 3.30, 3.65, 4.00,
        0.50, 0.85, 1.20, 1.55, 1.90, 2.25, 2.60, 2.95, 3.30, 3.65, 4.00,
        0.50, 0.85, 1.20, 1.55, 1.90, 2.25, 2.60, 2.95, 3.30, 3.65, 4.00,
        0.50, 0.85, 1.20, 1.55, 1.90, 2.25, 2.60, 2.95, 3.30, 3.65, 4.00,
        0.50, 0.85, 1.20, 1.55, 1.90, 2.25, 2.60, 2.95, 3.30, 3.65, 4.00,
        0.50, 0.85, 1.20, 1.55, 1.90, 2.25, 2.60, 2.95, 3.30, 3.65, 4.00,
        0.50, 0.85, 1.20, 1.55, 1.90, 2.25, 2.60, 2.95, 3.30, 3.65, 4.00,
        0.50, 0.85, 1.20, 1.55, 1.90, 2.25, 2.60, 2.95, 3.30, 3.65, 4.00,
        0.50, 0.85, 1.20, 1.55, 1.90, 2.25, 2.60, 2.95, 3.30, 3.65, 4.00,
        0.50, 0.85, 1.20, 1.55, 1.90, 2.25, 2.60, 2.95, 3.30, 3.65, 4.00,
        0.50, 0.85, 1.20, 1.55, 1.90, 2.25, 2.60, 2.95, 3.30, 3.65, 4.00,
	};

	PyObject* args = PyTuple_New(4);
	// c
	PyTuple_SET_ITEM(args, 0, PyArray_New(&PyArray_Type, 2, new npy_intp[2]{ 11, 11 }, NPY_DOUBLE, NULL, (double*)& c_, 0, NPY_ARRAY_CARRAY, NULL));
	// x0
	PyTuple_SET_ITEM(args, 1, PyArray_New(&PyArray_Type, 1, new npy_intp[1]{ 2 }, NPY_LONG, NULL, new long[2]{ 2, 0 }, 0, NPY_ARRAY_CARRAY, NULL));
	// dx
	PyTuple_SET_ITEM(args, 2, PyArray_New(&PyArray_Type, 1, new npy_intp[1]{ 2 }, NPY_DOUBLE, NULL, new double[2]{ 0.4, 0.7 }, 0, NPY_ARRAY_CARRAY, NULL));
	// order
	PyTuple_SET_ITEM(args, 3, PyLong_FromLong(2));

	PyArrayObject* tau = (PyArrayObject*) fast_marching(NULL, args);
	double* tau_c = (double*)PyArray_DATA(tau);

	for (int i = 0; i < 121; i++)
	{
		if (i > 0 && i % 11 == 0)
			cout << endl;
		printf("%2.2f\t", tau_c[i]);
	}

	Py_Finalize();
	//_CrtDumpMemoryLeaks();

	return 0;
}