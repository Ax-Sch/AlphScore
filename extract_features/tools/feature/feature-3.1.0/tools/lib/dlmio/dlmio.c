// $Id: dlmio.c,v 1.1.1.1 2004/05/22 01:18:15 mliang Exp $
// Copyright (c) 2004 Mike Liang. All rights reserved.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#undef DEBUG

#define BUF_SIZE 65536

double *dlmreadf(FILE *file,double **_values,int *_numrows,int *_numcols,char *delim,int skiprows,int skipcols) {
    long startpos;
    char buf[BUF_SIZE];
    char *token;
    double value;
    double *values;
    double *cursor;
    int numrows=0,numcols=0;
    int row=0,col=0;
    int skiprowidx, skipcolidx;


    /* count size */
    startpos = ftell(file);
    numrows = 0;
    // skip rows
    for (skiprowidx = 0; skiprowidx < skiprows; skiprowidx++) {
        if (fgets(buf,BUF_SIZE,file) == NULL) {
            break;
        }
    }
    // count cols
    if (fgets(buf,BUF_SIZE,file) != NULL) {
        numcols = 0;
        skipcolidx = 0;
        token = buf;
        do {
            // skip colS
            if (skipcolidx < skipcols) {
                skipcolidx++;
            } else {
                numcols++;
            }
            token = strstr(token+1,delim);
        } while (token);
        numrows++;
    }
    // count rows
    while (fgets(buf,BUF_SIZE,file) != NULL) {
        numrows++;
    }

#ifdef DEBUG
    fprintf(stderr,"Rows:%d Cols:%d\n",numrows,numcols);
#endif
    values = (double *)malloc(numrows*numcols*sizeof(double));
    if (!values) {
        fprintf(stderr,"Error allocating memory\n");
        return NULL;
    }

    /* read values */
    fseek(file,startpos,SEEK_SET);
    cursor = values;
    row = 0;
    for (skiprowidx = 0; skiprowidx < skiprows; skiprowidx++) {
        if (fgets(buf,BUF_SIZE,file) == NULL) {
            break;
        }
    }
    while (fgets(buf,BUF_SIZE,file) != NULL) {
        col = 0;
        skipcolidx = 0;
        token = buf;
        do {
            if (skipcolidx < skipcols) {
                skipcolidx++;
            } else {
                value = atof(token);
                *cursor++ = value;
                col++;
            }
            token = strstr(token+1,delim);
        } while (token);
        if (col != numcols) {
            fprintf(stderr,"Cols don't match %d != %d\n",col,numcols);
        }
        row++;
    }
    if (row != numrows) {
        fprintf(stderr,"Rows don't match %d != %d\n",row,numrows);
    }

    *_values = values;
    *_numrows = numrows;
    *_numcols = numcols;
    return values;
}

double *dlmread(char *filename,double **_values,int *_numrows,int *_numcols,char *delim,int skiprows,int skipcols) {
    double *retval;
    FILE *file;
    file = fopen(filename,"r");
    if (!file) {
        fprintf(stderr,"Error opening file %s\n",filename);
        return NULL;
    }
    retval = dlmreadf(file,_values,_numrows,_numcols,delim,skiprows,skipcols);
    fclose(file);
    return retval;
}

// ===================================================================
// PYTHON BINDING
// ===================================================================

#include "Python.h"
#include "Numeric/arrayobject.h"

static PyObject *
dlmio_dlmread(PyObject *self, PyObject *args, PyObject *kwargs) {
    static char *kwlist[] = {"fileobj","delim","row","col",NULL};
    PyObject *fileobj;
    FILE *file;
    char *filename;
    char *delim = "\t";
    double *values;
    PyArrayObject *array;
    int dims[2]; /* { rows, cols } */
    int skiprows=0;
    int skipcols=0;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|sii:dlmread", kwlist,
        &fileobj, &delim, &skiprows, &skipcols))
        return NULL;
    
    if (PyString_Check(fileobj)) {
        filename = PyString_AsString(fileobj);
        if (dlmread(filename,&values,&dims[0],&dims[1],delim,skiprows,skipcols) == NULL) {
            PyErr_SetFromErrno(PyExc_IOError);
            return NULL;
        }
    } else if (PyFile_Check(fileobj)) {
        file = PyFile_AsFile(fileobj);
        if (dlmreadf(file,&values,&dims[0],&dims[1],delim,skiprows,skipcols) == NULL) {
            PyErr_SetFromErrno(PyExc_IOError);
            return NULL;
        }
    } else {
        PyErr_SetString(PyExc_ValueError,"invalid file type");
        return NULL;
    }

#ifdef DEBUG
    fprintf(stderr,"Creating PyArray\n");
#endif
    array = (PyArrayObject *)PyArray_FromDims(2,dims,PyArray_DOUBLE);
    memcpy(array->data,values,dims[0]*dims[1]*sizeof(double));
    free(values);
    
    return PyArray_Return(array);
}

static PyMethodDef dlmioMethods[] = {
    {"dlmread", (PyCFunction)dlmio_dlmread, METH_VARARGS|METH_KEYWORDS, 
     "Delimited read\n"
     "\n"
     "dlmread(filename,delim='\\t',row=0,col=0)\n"
     "\n"
     "    filename - string or file object\n"
     "    delim - field delimiter\n"
     "    row - number of rows to skip\n"
     "    col - number of columns to skip\n"},
    {NULL,NULL,0,NULL} /* Sentinel */
};

PyMODINIT_FUNC
initdlmio(void) {
    (void) Py_InitModule("dlmio", dlmioMethods);
    import_array();
}
