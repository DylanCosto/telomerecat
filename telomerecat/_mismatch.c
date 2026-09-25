/* Compare sequence offsets using the public CPython API only. */
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdint.h>

static PyObject *compare_to_telo(PyObject *self, PyObject *args)
{
    PyObject *sequence, *quality, *pattern, *best = NULL;
    const unsigned char *seq, *qual, *pat;
    Py_ssize_t length, pattern_length, offset, index, position, score, best_score;
    uint64_t quality_sum;
    double average, best_average = Py_HUGE_VAL;
    unsigned int signal_counter = 0;

    if (!PyArg_ParseTuple(args, "OOO:compare_to_telo", &sequence, &quality, &pattern))
        return NULL;
    /* Unusual inputs retain the Python implementation and its exceptions. */
    if (!PyUnicode_CheckExact(sequence) || !PyUnicode_CheckExact(quality) ||
        !PyUnicode_CheckExact(pattern) || !PyUnicode_IS_ASCII(sequence) ||
        !PyUnicode_IS_ASCII(quality) || !PyUnicode_IS_ASCII(pattern))
        Py_RETURN_NOTIMPLEMENTED;
    length = PyUnicode_GET_LENGTH(sequence);
    pattern_length = PyUnicode_GET_LENGTH(pattern);
    if (PyUnicode_GET_LENGTH(quality) != length || pattern_length == 0 ||
        (uint64_t)length > INT64_MAX / 127)
        Py_RETURN_NOTIMPLEMENTED;
    seq = PyUnicode_1BYTE_DATA(sequence);
    qual = PyUnicode_1BYTE_DATA(quality);
    pat = PyUnicode_1BYTE_DATA(pattern);
    best_score = length;

    for (offset = 0; offset < pattern_length; ++offset) {
        score = 0;
        quality_sum = 0;
        position = offset;
        for (index = 0; index < length; ++index) {
            if (++signal_counter == 65536) {
                signal_counter = 0;
                if (PyErr_CheckSignals() < 0)
                    goto error;
            }
            if (seq[index] != pat[position]) {
                quality_sum += qual[index];
                if (++score > best_score)
                    break;
            }
            if (++position == pattern_length)
                position = 0;
        }
        if (score > best_score)
            continue;
        average = score ? (double)quality_sum / (double)score : 0.0;
        /* Strict comparison preserves the first offset on an exact tie. */
        if (best && score == best_score && average >= best_average)
            continue;
        PyObject *loci = PyList_New(score);
        if (!loci)
            goto error;
        position = offset;
        Py_ssize_t next = 0;
        for (index = 0; index < length; ++index) {
            if (seq[index] != pat[position]) {
                PyObject *locus = PyLong_FromSsize_t(index);
                if (!locus) {
                    Py_DECREF(loci);
                    goto error;
                }
                PyList_SET_ITEM(loci, next++, locus);
            }
            if (++position == pattern_length)
                position = 0;
        }
        Py_XDECREF(best);
        best = loci;
        best_score = score;
        best_average = average;
    }
    return best;

error:
    Py_XDECREF(best);
    return NULL;
}

static PyMethodDef methods[] = {
    {"compare_to_telo", compare_to_telo, METH_VARARGS,
     "Return mismatch loci for the lowest-count, lowest-mean-quality offset."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef module = {
    PyModuleDef_HEAD_INIT, "_mismatch", NULL, -1, methods
};

PyMODINIT_FUNC PyInit__mismatch(void)
{
    return PyModule_Create(&module);
}
