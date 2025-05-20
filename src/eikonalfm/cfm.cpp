// interface between Python and the C++ implementation
// https://docs.microsoft.com/de-de/visualstudio/python/working-with-c-cpp-python-in-visual-studio?view=vs-2019
// https://docs.scipy.org/doc/numpy-1.13.0/user/c-info.how-to-extend.html

#ifndef Py_cfm
#define Py_cfm

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

// Make "s#" use Py_ssize_t rather than int.
#define PY_SSIZE_T_CLEAN
#ifdef _DEBUG
	#undef _DEBUG
	#include <Python.h>
	#define _DEBUG
#else
	#include <Python.h>
#endif

#include <stdexcept>
#include <cstdio>
#include <signal.h>
#include <memory>
#include <vector>
// #include "numpy/noprefix.h" includes "no-prefix option" of the types and methods, i.e. intp instead of npy_intp
// defining NPY_NO_PREFIX prior to including "numpy/arrayobject.h" does this as well
#include "numpy/arrayobject.h"

#include "factoredmarcher.hpp"

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

        objects.push_back((PyObject *)obj);
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
    int size_s = std::snprintf(nullptr, 0, format_str, args...);
    if(size_s <= 0){ throw std::runtime_error("Error during formatting."); }
    auto size = static_cast<size_t>(size_s);
    std::unique_ptr<char[]> msg(new char[size]);
    std::snprintf(msg.get(), size_s, format_str, args...);
    return msg;
}

// we don't want the function names to get 'mangled' by the compiler
#ifdef __cplusplus
extern "C" {
#endif

static PyObject *fast_marching_(PyObject *args, PyObject *kwargs, const bool factored)
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
	usize *shape = new usize[ndim];
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

	delete[] shape;
	delete m;

    // restore original signal handler
    signal(SIGINT, handler_sigint);

	if (!success)
	{
        collection.DECREF();
        return NULL;
	}

	
	PyObject *res;
	if (output_sensitivities)
	{
		npy_intp *npy_shape = PyArray_DIMS(c);

		npy_intp *orders_shape = new npy_intp[ndim + 1];
		orders_shape[0] = ndim; // one array for each dimension
		for (int d=0; d < ndim; d++)
			orders_shape[d+1] = npy_shape[d];

		// the pointers of `info` must not be deallocated! otherwise the numpy array is broken
		PyObject *sequence = PyArray_SimpleNewFromData(ndim, npy_shape, NPY_LONG, info->get_sequence());
		if (!collection.validate(sequence, PyExc_MemoryError, "couldn't create sequence array."))
		{
			delete[] info->get_sequence();
			delete[] info->get_order();
			return NULL;
		}
		// set the array object to own the data pointer and therefore call delete[] when it is destroyed
		PyArray_ENABLEFLAGS((PyArrayObject *)sequence, NPY_ARRAY_OWNDATA);

		PyObject *orders = PyArray_SimpleNewFromData(ndim+1, orders_shape, NPY_INT8, info->get_order());
		if (!collection.validate(orders, PyExc_MemoryError, "couldn't create sequence array."))
		{
			// info->get_sequence is already owned by the array, so don't delete it here
			delete[] info->get_order();
			return NULL;
		}
		PyArray_ENABLEFLAGS((PyArrayObject *)orders, NPY_ARRAY_OWNDATA);

		// this will increase the ref counter of each object by 1!
		res = Py_BuildValue("OOO", tau, sequence, orders);

		delete[] orders_shape;
	}
	else
	{
		Py_INCREF(tau); // need to increase by 1 here so that after collection.DECREF the count is still 1
		res = (PyObject *)tau;
	}

	collection.DECREF();
	return res;
}

static PyObject *fast_marching_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
	return fast_marching_(args, kwargs, false);
}

static PyObject *factored_marching_wrapper(PyObject *self, PyObject *args, PyObject *kwargs)
{
	return fast_marching_(args, kwargs, true);
}

// Each entry in this array is a PyMethodDef structure containing 
// 1) the Python name,
// 2) the C - function that implements the function,
// 3) flags indicating whether or not keywords are accepted for this function, and 
// 4) The docstring for the function.
static PyMethodDef fm_methods[] = {
	{"fast_marching", (PyCFunction)fast_marching_wrapper, METH_VARARGS | METH_KEYWORDS,
        "fast_marching(c, x_s, dx, order, output_sensitivities=False)\n--\n\n"
        "    Calculates the fast marching solution to the eikonal equation.\n"
        "\n"
        "    Parameters\n"
        "    ----------\n"
        "    c : ndarray\n"
        "        (background) velocity array, c(x) > 0.\n"
        "    x_s : sequence of ints\n"
        "        Source position as index vector, e.g. ``(0, 0)``. Must have the same length as the number of dimensions of c.\n"
        "    dx : sequence of doubles\n"
        "        Grid spacing for each dimension, dx > 0. Must have the same length as the number of dimensions of c.\n"
        "    order : {1, 2}\n"
        "        Order of the finite difference operators.\n"
		"    output_sensitivities : boolean, optional\n"
		"        Additionally returns sensitivity data if True.\n"
        "\n"
        "    Returns\n"
        "    ----------\n"
        "    tau : ndarray\n"
        "        numerical solution tau for the eikonal equation.\n"
		"    sequence : ndarray, optional\n"
		"        Only returned when 'output_sensitivities=True'. The sequence in which each gridpoint was set to known.\n"
		"    orders : ndarray, optional\n"
		"        Only returned when 'output_sensitivities=True'. The finite difference orders for each dimension used at each gridpoint."
    },
	{"factored_fast_marching", (PyCFunction)factored_marching_wrapper, METH_VARARGS | METH_KEYWORDS,
        "factored_fast_marching(c, x_s, dx, order, output_sensitivities=False)\n--\n\n"
        "    Calculates the fast marching solution to the factored eikonal equation.\n"
        "\n"
        "    Parameters\n"
        "    ----------\n"
        "    c : ndarray\n"
        "        (background) velocity array, c(x) > 0.\n"
        "    x_s : sequence of ints\n"
        "        Source position as index vector, e.g. ``(0, 0)``. Must have the same length as the number of dimensions of c.\n"
        "    dx : sequence of doubles\n"
        "        Grid spacing for each dimension, dx > 0. Must have the same length as the number of dimensions of c.\n"
        "    order : {1, 2}\n"
        "        Order of the finite difference operators.\n"
		"    output_sensitivities : boolean, optional\n"
		"        Additionally returns sensitivity data if True.\n"
        "\n"
        "    Returns\n"
        "    ----------\n"
        "    tau1 : ndarray\n"
        "        numerical solution tau1 for the factored eikonal equation.\n"
        "        To get tau, you need to multiply this with the distance field tau0.\n"
		"    sequence : ndarray, optional\n"
		"        Only returned when 'output_sensitivities=True'. The sequence in which each gridpoint was set to known.\n"
		"    orders : ndarray, optional\n"
		"        Only returned when 'output_sensitivities=True'. The finite difference orders for each dimension used at each gridpoint."
    },

	// Terminate the array with an object containing nulls.
	{NULL, NULL, 0, NULL}
};

static struct PyModuleDef cfm_module = {
	PyModuleDef_HEAD_INIT,
	"cfm",                                                                 		// Module name to use with Python import statements
	"c++ implementation of (factored) fast marching for the eikonal equation",  // Module description
	0,																			// size of per-interpreter state of the module, or -1 if the module keeps state in global variables.
	fm_methods																	// Structure that defines the methods of the module
};

PyMODINIT_FUNC PyInit_cfm(void)
{
	PyObject *m = PyModule_Create(&cfm_module);
	if (m == NULL)
		return NULL;

	import_array();
	return m;
}

#ifdef __cplusplus
}
#endif

#endif // !defined(Py_cfm)