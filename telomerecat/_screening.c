/* Screen paired reads using pysam's public Python interface.
 * No pysam or HTSlib headers, private structures, or binary ABI are needed.
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <string.h>

typedef struct {
    PyObject *first;
    PyObject *second;
    const char *first_text;
    const char *second_text;
    int use_strings;
} Patterns;

static int prepare_patterns(Patterns *patterns)
{
    Py_ssize_t first_length, second_length;
    patterns->use_strings = 0;
    if (!PyUnicode_Check(patterns->first) || !PyUnicode_Check(patterns->second) ||
        !PyUnicode_IS_ASCII(patterns->first) || !PyUnicode_IS_ASCII(patterns->second))
        return 0;
    patterns->first_text = PyUnicode_AsUTF8AndSize(patterns->first, &first_length);
    if (!patterns->first_text)
        return -1;
    patterns->second_text = PyUnicode_AsUTF8AndSize(patterns->second, &second_length);
    if (!patterns->second_text)
        return -1;
    patterns->use_strings =
        memchr(patterns->first_text, 0, first_length) == NULL &&
        memchr(patterns->second_text, 0, second_length) == NULL;
    return 0;
}

static int matches_read(PyObject *read, Patterns *patterns)
{
    PyObject *sequence = PyObject_GetAttrString(read, "query_sequence");
    Py_ssize_t length;
    const char *text;
    int found;
    if (!sequence)
        return -1;
    if (patterns->use_strings && PyUnicode_Check(sequence) && PyUnicode_IS_ASCII(sequence)) {
        text = PyUnicode_AsUTF8AndSize(sequence, &length);
        if (!text) {
            Py_DECREF(sequence);
            return -1;
        }
        /* BAM sequences cannot contain NUL, but preserve Python containment
         * semantics for other objects exposing this same interface too.
         */
        if (memchr(text, 0, length) == NULL) {
            found = strstr(text, patterns->first_text) != NULL ||
                    strstr(text, patterns->second_text) != NULL;
            Py_DECREF(sequence);
            return found;
        }
    }
    found = PySequence_Contains(sequence, patterns->first);
    if (found == 0)
        found = PySequence_Contains(sequence, patterns->second);
    Py_DECREF(sequence);
    return found;
}

static PyObject *pairs_to_telbam(PyObject *self, PyObject *args)
{
    PyObject *source, *destination;
    PyObject *fetch = NULL, *empty_args = NULL, *kwargs = NULL;
    PyObject *fetched = NULL, *iterator = NULL, *write = NULL;
    PyObject *first = NULL, *second = NULL, *written = NULL;
    Patterns patterns;
    unsigned int since_signal_check = 0;
    int found;
    if (!PyArg_ParseTuple(args, "OOOO:pairs_to_telbam", &source, &destination,
                          &patterns.first, &patterns.second))
        return NULL;
    if (prepare_patterns(&patterns) < 0)
        return NULL;
    fetch = PyObject_GetAttrString(source, "fetch");
    if (!fetch)
        goto error;
    empty_args = PyTuple_New(0);
    if (!empty_args)
        goto error;
    kwargs = Py_BuildValue("{s:O}", "until_eof", Py_True);
    if (!kwargs)
        goto error;
    fetched = PyObject_Call(fetch, empty_args, kwargs);
    if (!fetched)
        goto error;
    iterator = PyObject_GetIter(fetched);
    if (!iterator)
        goto error;
    Py_CLEAR(fetch);
    Py_CLEAR(empty_args);
    Py_CLEAR(kwargs);
    Py_CLEAR(fetched);
    write = PyObject_GetAttrString(destination, "write");
    if (!write)
        goto error;

    while (1) {
        first = PyIter_Next(iterator);
        if (!first) {
            if (PyErr_Occurred())
                goto error;
            break;
        }
        second = PyIter_Next(iterator);
        if (!second) {
            if (!PyErr_Occurred())
                PyErr_SetNone(PyExc_StopIteration);
            goto error;
        }
        found = matches_read(first, &patterns);
        if (found == 0)
            found = matches_read(second, &patterns);
        if (found < 0)
            goto error;
        if (found) {
            written = PyObject_CallFunctionObjArgs(write, first, NULL);
            if (!written)
                goto error;
            Py_CLEAR(written);
            written = PyObject_CallFunctionObjArgs(write, second, NULL);
            if (!written)
                goto error;
            Py_CLEAR(written);
        }
        Py_CLEAR(first);
        Py_CLEAR(second);
        if (++since_signal_check == 65536) {
            if (PyErr_CheckSignals() < 0)
                goto error;
            since_signal_check = 0;
        }
    }
    Py_DECREF(iterator);
    Py_DECREF(write);
    Py_RETURN_NONE;

error:
    Py_XDECREF(fetch);
    Py_XDECREF(empty_args);
    Py_XDECREF(kwargs);
    Py_XDECREF(fetched);
    Py_XDECREF(iterator);
    Py_XDECREF(write);
    Py_XDECREF(first);
    Py_XDECREF(second);
    Py_XDECREF(written);
    return NULL;
}

static PyMethodDef methods[] = {
    {"pairs_to_telbam", pairs_to_telbam, METH_VARARGS,
     "Write both reads, in input order, when either contains either motif."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef module = {
    PyModuleDef_HEAD_INIT, "_screening", NULL, -1, methods
};

PyMODINIT_FUNC PyInit__screening(void)
{
    return PyModule_Create(&module);
}
