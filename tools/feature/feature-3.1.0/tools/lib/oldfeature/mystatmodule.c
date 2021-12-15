/* $Id: mystatmodule.c,v 1.1 2001/12/04 19:23:25 mliang Exp $ */
/* Copyright (c) 2001 Mike Liang.  All rights reserved. */

/* Python wrapper to statistic functions in C */

#include "Python.h"
#include "normal_cdf.h"



/************************************** Exported Functions ***********/
static char mystat_normal_cdf__doc__[] = 
"normal_cdf(z)";

static PyObject *
mystat_normal_cdf(self, args)
     PyObject *self;
     PyObject *args;
{
    double z, p;
    if(!PyArg_ParseTuple(args, "d", &z))
	return NULL;
    p = normal_cdf(z);
    return PyFloat_FromDouble(p);
}


/************************************** Module definition stuff ******/

static PyMethodDef mystatMethods[] = {
    {"normal_cdf", mystat_normal_cdf, METH_VARARGS, mystat_normal_cdf__doc__},
    {NULL, NULL}
};

static char mystat__doc__[] =
"My Statistic functions implemented in C.\n";

void initmystat()
{
    (void) Py_InitModule3("mystat", mystatMethods, mystat__doc__);
}
