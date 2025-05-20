https://stackoverflow.com/questions/52731884/pyarray-simplenewfromdata

The technique is a good one: create a Python object that knows how to free your memory when destroyed, and make it the base of the returned array.

Python comes with a type designed to perform arbitrary C-level cleanup when destroyed: capsules, which bundle together a pointer and a destructor function and call the destructor when the capsule is destroyed.

To create a capsule for your memory, first, we define a destructor function:

```c++
void capsule_cleanup(PyObject *capsule) {
    void *memory = PyCapsule_GetPointer(capsule, NULL);
    // I'm going to assume your memory needs to be freed with free().
    // If it needs different cleanup, perform whatever that cleanup is
    // instead of calling free().
    free(memory);
}
```

And you set a capsule as your array's base with

```c++
PyObject *capsule = PyCapsule_New(data, NULL, capsule_cleanup);
PyArray_SetBaseObject((PyArrayObject *) arr, capsule);
// Do not Py_DECREF the capsule; PyArray_SetBaseObject stole your
// reference.
```

And that should ensure your memory gets freed once it's no longer in use.
