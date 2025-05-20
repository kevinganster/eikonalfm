// memory leak detection
//#include "debugtests.h"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#define PY_SSIZE_T_CLEAN
#ifdef _DEBUG
	#undef _DEBUG
	#include <Python.h>
	#define _DEBUG
#else
	#include <Python.h>
#endif

#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cstdio>
#include <signal.h>
#include <memory>
#include <vector>
#include "numpy/arrayobject.h"
#include "factoredmarcher.hpp"


// copied from cfm.cpp
class InterruptExc : public std::exception
{
public:
    InterruptExc() {}
};

void signal_handler(int code)
{
   throw InterruptExc();
};

class PyObjCollection
{
private:
    std::vector<PyObject *> objects;

public:
    /**
     * Validates an PyObject or PyArrayObject.
     * If it fails, a specified error will be returned to Python and all stored pointers will be dereferenced.
     * Otherwise it adds the pointer to the list of stored pointers.
     */
    template <class T>
    bool validate(T *obj, PyObject *error, const char *msg)
    {
        if (!obj)
        {
            PyErr_SetString(error, msg);
            DECREF();
            return false;
        }

        objects.push_back(_PyObject_CAST(obj));
        return true;
    }

    void DECREF()
    {
        for (auto py_obj : objects)
        {
            Py_DECREF(py_obj);
        }
    }
};

template <class... Args>
auto format(const char *format_str, Args&&... args) -> std::unique_ptr<char[]>
{
    int size_s = std::snprintf(nullptr, 0, format_str, args ...);
    if(size_s <= 0){ throw std::runtime_error("Error during formatting."); }
    auto size = static_cast<size_t>(size_s);
    std::unique_ptr<char[]> msg(new char[size]);
    std::sprintf(msg.get(), format_str, args...);
    return msg;
}

static PyObject *fast_marching(PyObject *args, PyObject *kwargs, const bool factored)
{
    // will be const_cast<char**>() below, hopefully PyArg_ParseTupleAndKeywords doesn't try to alter it
    const char *keywords[] = {"c", "x_s", "dx", "order", "output_sensitivities", NULL};
	// define placeholders
	PyObject *pc, *px_s, *pdx;
	int order;
	bool output_sensitivities = false;

	// parse arguments
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOi|b", const_cast<char**>(keywords), &pc, &px_s, &pdx, &order, &output_sensitivities))
		return NULL;

	// only orders 1 and 2 are implemented
	if (order < 1 || order > 2)
	{
		PyErr_SetString(PyExc_ValueError, "Only orders 1 and 2 are supported.");
		return NULL;
	}

	PyObjCollection collection;

	// convert python-numpy-array to c-numpy-array
	PyArrayObject *c = (PyArrayObject *)PyArray_FROMANY(pc, NPY_DOUBLE, 1, 0, NPY_ARRAY_IN_ARRAY);
	if (!collection.validate(c, PyExc_ValueError, "c must be an array of doubles."))
		return NULL;

	int ndim = PyArray_NDIM(c);

	// check if px0 is a python-numpy-array and if it has the correct size
	PyArrayObject *x_s_ = (PyArrayObject *)PyArray_FROMANY(px_s, NPY_LONG, 1, 1, NPY_ARRAY_IN_ARRAY);
	if (!collection.validate(x_s_, PyExc_ValueError, "x_s must be a 1D array of nonnegative ints."))
		return NULL;
	else if (PyArray_SIZE(x_s_) != ndim)
	{
		PyErr_SetString(PyExc_ValueError, format("size of x_s and number of dimensions of c do not match: %d != %d.", (int)PyArray_SIZE(x_s_), ndim).get());
        collection.DECREF();
		return NULL;
	}

	// check if pdx is a python-numpy-array and if it has the correct size
	PyArrayObject *dx_ = (PyArrayObject *)PyArray_FROMANY(pdx, NPY_DOUBLE, 1, 1, NPY_ARRAY_IN_ARRAY);
	if (!collection.validate(dx_, PyExc_ValueError, "dx must be a 1D array of doubles."))
		return NULL;
	else if (PyArray_SIZE(dx_) != ndim)
	{
		PyErr_SetString(PyExc_ValueError, format("size of dx and number of dimensions of c do not match: %d != %d.", (int)PyArray_SIZE(dx_), ndim).get());
        collection.DECREF();
		return NULL;
	}

	double *dx = (double*)PyArray_DATA(dx_);
	// TODO: maybe change this to 'ssize', since x_s_d is of type 'ssize' anyways
	usize shape[ndim];
	// index version of x_s
    usize x_s = 0;
    usize *x_s_d = (usize *)PyArray_DATA(x_s_);

    usize tmp = 1;
	for (int d = ndim - 1; d >= 0; d--)
	{
		shape[d] = PyArray_DIM(c, d);

		if (x_s_d[d] >= shape[d])
		{
			PyErr_SetString(PyExc_ValueError, format("entries of x_s have to be within the shape of c, but entry %d with value %ld is not in [0, %ld].", d, x_s_d[d], shape[d]-1).get());
			collection.DECREF();
			return NULL;
		}
		if (dx[d] <= 0)
		{
			PyErr_SetString(PyExc_ValueError, format("entries of dx must be greater than zero, but entry %d with value %f is not.", d, dx[d]).get());
			collection.DECREF();
			return NULL;
		}

		x_s += tmp * x_s_d[d];
		tmp *= shape[d];
	}

	// create a new array for the return value
	PyArrayObject *tau = (PyArrayObject *)PyArray_ZEROS(ndim, PyArray_DIMS(c), NPY_DOUBLE, 0);
	if (!collection.validate(tau, PyExc_MemoryError, "couldn't create result array."))
		return NULL;

	std::unique_ptr<MarcherInfo> info;
	if (output_sensitivities)
		info = std::unique_ptr<MarcherInfo>(new SensitivityInfo{ndim, shape});
	else
		info = std::unique_ptr<MarcherInfo>(new MarcherInfo{ndim, shape});


	Marcher *m;
	if (factored)
		m = new FactoredMarcher{(double *)PyArray_DATA(c), *info, dx, order};
	else
		m = new Marcher{(double *)PyArray_DATA(c), *info, dx, order};

    // install custom sigintterrupt handler (CTRL+C)
    auto handler_sigint = signal(SIGINT, signal_handler);
    bool success = true;

	try
	{
		m->solve(x_s, (double *)PyArray_DATA(tau));
	}
	catch (const InterruptExc& ex)
    {
        PyErr_SetString(PyExc_KeyboardInterrupt, "");
        success = false;
    }
	catch (const std::exception& ex)
    {
        // propagate error
        PyErr_SetString(PyExc_RuntimeError, ex.what());
        success = false;
    }

	delete m;

    // restore original signal handler
    signal(SIGINT, handler_sigint);

	if (!success)
	{
        collection.DECREF();
        return NULL;
	}

	
	if (output_sensitivities)
	{
		npy_intp *npy_shape = PyArray_DIMS(c);

		npy_intp orders_shape[ndim + 1];
		orders_shape[0] = ndim; // one array for each dimension
		for (int d=0; d < ndim; d++)
			orders_shape[d+1] = npy_shape[d];

		PyObject *sequence = PyArray_New(&PyArray_Type, ndim, npy_shape, NPY_LONG, NULL, info->get_sequence(), 0, NPY_ARRAY_CARRAY, NULL);
		if (!collection.validate(sequence, PyExc_MemoryError, "couldn't create sequence array."))
			return NULL;
		PyObject *orders = PyArray_New(&PyArray_Type, ndim+1, orders_shape, NPY_INT8, NULL, info->get_order(), 0, NPY_ARRAY_CARRAY, NULL);
		if (!collection.validate(orders, PyExc_MemoryError, "couldn't create sequence array."))
			return NULL;
		return Py_BuildValue("OOO", (PyObject *)tau, sequence, orders);
	}
	else
		return (PyObject *)tau;
}


int main()
{
	//Py_SetPath(Py_DecodeLocale("E:\\Anaconda\\Lib;E:\\Anaconda\\Lib\\site-packages;E:\\Anaconda\\DLLs", NULL));
    Py_SetPath(Py_DecodeLocale("~\\anaconda3\\lib;~\\anaconda3\\lib\\python3.9\\site-packages;", NULL));
    // Py_SetPath(Py_DecodeLocale("D:\\Funktionsprogramme\\Anaconda\\lib;D:\\Funktionsprogramme\\Anaconda\\lib\\site-packages;D:\\Funktionsprogramme\\Anaconda\\DLLs", NULL));

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

	PyObject *args = PyTuple_New(4);
	// c
	PyTuple_SET_ITEM(args, 0, PyArray_New(&PyArray_Type, 2, new npy_intp[2]{ 11, 11 }, NPY_DOUBLE, NULL, (double *)& c_, 0, NPY_ARRAY_CARRAY, NULL));
	// x_s
	PyTuple_SET_ITEM(args, 1, PyArray_New(&PyArray_Type, 1, new npy_intp[1]{ 2 }, NPY_LONG, NULL, new long[2]{ 2, 0 }, 0, NPY_ARRAY_CARRAY, NULL));
	// dx
	PyTuple_SET_ITEM(args, 2, PyArray_New(&PyArray_Type, 1, new npy_intp[1]{ 2 }, NPY_DOUBLE, NULL, new double[2]{ 0.4, 0.7 }, 0, NPY_ARRAY_CARRAY, NULL));
	// order
	PyTuple_SET_ITEM(args, 3, PyLong_FromLong(2));

	PyObject *kwargs = PyDict_New();

	PyArrayObject *tau = (PyArrayObject *) fast_marching(args, kwargs, false);
	double *tau_c = (double *)PyArray_DATA(tau);

	for (int i = 0; i < 121; i++)
	{
		if (i > 0 && i % 11 == 0)
			std::cout << std::endl;
		printf("%2.2f\t", tau_c[i]);
	}

	Py_Finalize();
	//_CrtDumpMemoryLeaks();

	return 0;
}