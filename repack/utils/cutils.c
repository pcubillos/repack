// Copyright (c) 2017-2020 Patricio Cubillos and contributors.
// repack is open-source software under the MIT license (see LICENSE).

#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <math.h>

/* 1D double ndarray:                                                       */
#define INDd(a,i) *((double *)(PyArray_DATA(a) + i * PyArray_STRIDE(a, 0)))
/* 1D integer ndarray:                                                      */
#define INDi(a,i) *((int    *)(PyArray_DATA(a) + i * PyArray_STRIDE(a, 0)))
#define INDb(a,i) *((_Bool  *)(PyArray_DATA(a) + i * PyArray_STRIDE(a, 0)))


int
bsfind(PyArrayObject *array, double value, int lo, int hi){
  /* Last case, value limited between consecutive indices of array:         */
  if (hi-lo <= 1){
    /* Return closest array index to value:                                 */
    if (fabs(INDd(array,hi)-value) < fabs(INDd(array,lo)-value))
      return hi;
    return lo;
  }
  /* Compare to middle point and search in corresponding sub-array:         */
  else if (INDd(array,((hi+lo)/2)) > value)
    return bsfind(array, value, lo, (hi+lo)/2);
  else
    return bsfind(array, value, (hi+lo)/2, hi);
}


PyDoc_STRVAR(continuum__doc__,
"Compute tabulated extinction coefficient (cm2 molec-1)         \n\
by diluting line-transitions into wavenumber array.             \n\
                                                                \n\
Parameters                                                      \n\
----------                                                      \n\
s: 1D float ndarray                                             \n\
   Line-transition strength per Doppler width.                  \n\
wn: 1D float ndarray                                            \n\
   Central wavenumber of the line transitions.                  \n\
cont: 1D float ndarray                                          \n\
   Output extinction-coefficient values corresponding to wnspec.\n\
wnspec: 1D float ndarray                                        \n\
   Tabulated wavenumber array.                                  \n\
");

static PyObject *continuum(PyObject *self, PyObject *args){
  PyArrayObject *wn, *s, *wnspec, *cont;
  int i, j,           /* Auxilliary for-loop indices                        */
      nlines, nwave;  /* Number of lines and wavenumber samples             */
  double c, dwn;

  /* Load inputs:                                                           */
  if (!PyArg_ParseTuple(args, "OOOO", &s, &wn, &cont, &wnspec))
    return NULL;

  /* Get the number of lines and spectrum size:                             */
  nlines = (int)PyArray_DIM(s,     0);
  nwave  = (int)PyArray_DIM(wnspec,0);
  /* Wavenumber sampling rate:                                              */
  dwn = INDd(wnspec,1) - INDd(wnspec,0);
  /* Co-add line strengths to tabulated nearest neighbor wavenumber sample: */
  j = 0;
  for (i=0; i<nwave; i++){
    c = 0.0;
    while (INDd(wn,j) < INDd(wnspec,i)+0.5*dwn){
      if (j == nlines)
        break;
      c += INDd(s,j);
      j++;
    }
    /* and dilute by the wavenumber bin width:                              */
    INDd(cont,i) += c/dwn;
  }
  return Py_BuildValue("i", 1);
}


PyDoc_STRVAR(flag__doc__,
"Flag weak from strong line-transtion lines.                \n\
                                                            \n\
Parameters                                                  \n\
----------                                                  \n\
flag: 1D bool ndarray                                       \n\
   Output indicating the strong (1) and weak (0) lines.     \n\
wn: 1D float ndarray                                        \n\
   Central wavenumber of the line transitions.              \n\
s: 1D float ndarray                                         \n\
   Line-transition strength per Doppler width.              \n\
isort: 1D integer ndarray                                   \n\
   Indices of line transitions sorted in decreasing s order.\n\
alphad: 1D float ndarray                                    \n\
   Doppler width of the line transitions.                   \n\
sthreash: Float                                             \n\
   Tolerance threshold to define weak from strong lines     \n\
");

static PyObject *flag(PyObject *self, PyObject *args){
  PyArrayObject *flag, *wn, *s, *isort, *alphad;
  int i, j, k,              /* Auxilliary for-loop indices                  */
      nlines, imin, imax;
  double f, sthresh, dop;

  /* Load inputs:                                                           */
  if (!PyArg_ParseTuple(args, "OOOOOd", &flag, &wn, &s, &isort,
                                        &alphad, &sthresh))
    return NULL;

  /* Get the spectrum size:                                                 */
  nlines = (int)PyArray_DIM(s, 0);
  /* Evaluate the Planck function:                                          */
  for (j=0; j<nlines; j++){
    i = INDi(isort,j);
    if (INDb(flag,i)){
      /* Find limits                                                        */
      imin = bsfind(wn, INDd(wn,i)-6*INDd(alphad,i), 0, i);
      imax = bsfind(wn, INDd(wn,i)+6*INDd(alphad,i), i, nlines-1);
      /* Evaluate s against doppler:                                        */
      f = 1.0 / pow(INDd(alphad,i),2);
      for (k=imin; k<imax; k++){
        if (INDb(flag,k) && k != i){
          dop = exp(-f*pow(INDd(wn,k)-INDd(wn,i),2));
          INDb(flag,k) = INDd(s,k) > sthresh * INDd(s,i) * dop;
        }
      }
    }
  }
  return Py_BuildValue("i", 1);
}


/* The module doc string    */
PyDoc_STRVAR(cutils__doc__, "Repack C-utility functions.");

/* A list of all the methods defined by this module.                        */
static PyMethodDef cutils_methods[] = {
    {"continuum", continuum, METH_VARARGS, continuum__doc__},
    {"flag",      flag,      METH_VARARGS, flag__doc__},
    {NULL,        NULL,      0,            NULL}                /* sentinel */
};


#if PY_MAJOR_VERSION >= 3
/* Module definition for Python 3.                                          */
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "cutils",
    cutils__doc__,
    -1,
    cutils_methods
};

/* When Python 3 imports a C module named 'X' it loads the module           */
/* then looks for a method named "PyInit_"+X and calls it.                  */
PyObject *PyInit_cutils (void) {
  PyObject *module = PyModule_Create(&moduledef);
  import_array();
  return module;
}

#else
/* When Python 2 imports a C module named 'X' it loads the module           */
/* then looks for a method named "init"+X and calls it.                     */
void initcutils(void){
  Py_InitModule3("cutils", cutils_methods, cutils__doc__);
  import_array();
}
#endif
