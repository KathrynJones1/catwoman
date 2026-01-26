/* catwoman: a batman extension to generate morning/evening terminator transit lightcurves
 * Copyright (C) 2019 Kathryn Jones & Néstor Espinoza
 *
 * This program incorporates a modified version of the batman package: fast computation of exoplanet transit light curves
 * Copyright (C) 2015 Laura Kreidberg
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include "numpy/arrayobject.h"
#include "common.h"
#include <stdbool.h>

#if defined (_OPENACC) && defined(__PGI)
#  include <accelmath.h>
#else
#  include <math.h>
#endif

#if defined (_OPENMP) && !defined(_OPENACC)
#  include <omp.h>
#endif

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

inline double ellpic_bulirsch(double n, double k);
inline double ellec(double k);
inline double ellk(double k);
//double intensity(double x, double* args);
static PyObject *_quadratic_ld(PyObject *self, PyObject *args);


inline double intensity(double x, double* args)
{       
        double c1=args[0], c2=args[1];  //norm=args[4]; not sure what this is

        //if(x > 0.99995) x = 0.99995;  //ask Nestor if this should be kept in**
	double mu = sqrt(1.-x*x);
        return (1. - c1*(1. - mu) - c2*(1. - mu)*(1. - mu))/((1. - c1/3 - c2/6)*M_PI);
}


static PyObject *_quadratic_ld(PyObject *self, PyObject *args)
{
/*    Input: *************************************

     ds       array of impact parameters in units of rs
     c1       linear    limb-darkening coefficient (gamma_1 in Mandel & Agol 2002)
     c2       quadratic limb-darkening coefficient (gamma_2)
     p        occulting star(I think this is meant to be planet) size in units of rs

     Output: ***********************************

     flux     fraction of flux at each ds for a limb-darkened source

     Limb darkening has the form:
     I(r) = [1 - c1 * (1 - sqrt(1 - (r/rs)^2)) - c2*(1 - sqrt(1 - (r/rs)^2))^2]/(1 - c1/3 - c2/6)/pi
*/
    const int nthreads;
    const double c1, c2, p;
    PyArrayObject *ds, *flux, *phi, *b, *x;
    npy_intp dims[1];
    bool twoc;
    double fac = 0.0001;
    double rp2;

    if(!PyArg_ParseTuple(args,"OddddiOOOdb", &ds, &p, &c1, &c2, &fac, &nthreads, &phi, &b, &x, &rp2, &twoc)) return NULL;

    dims[0] = PyArray_DIMS(ds)[0];
    flux = (PyArrayObject *) PyArray_SimpleNew(1, dims, PyArray_TYPE(ds));    //creates numpy array to store return flux values

    double *f_array = PyArray_DATA(flux);
    double *d_array = PyArray_DATA(ds);
    
    double *b_array = PyArray_DATA(b);
    double *phi_array = PyArray_DATA(phi);
    double *x_array = PyArray_DATA(x);

    /*
    NOTE:  the safest way to access numpy arrays is to use the PyArray_GETITEM and PyArray_SETITEM functions.
    Here we use a trick for faster access and more convenient access, where we set a pointer to the
    beginning of the array with the PyArray_DATA (e.g., f_array) and access elements with e.g., f_array[i].
    Success of this operation depends on the numpy array storing data in blocks equal in size to a C double.
    If you run into trouble along these lines, I recommend changing the array access to something like:
        d = PyFloat_AsDouble(PyArray_GETITEM(ds, PyArray_GetPtr(ds, &i)));
    where ds is a numpy array object.
    Laura Kreidberg 07/2015
    */

	double intensity_args[] = {c1, c2}; //added this
    	#pragma acc data copyin(intensity_args)
    	calc_limb_darkening(f_array, d_array, dims[0], p, fac, nthreads, intensity_args, phi_array, b_array, x_array, rp2, twoc);
    
    	return PyArray_Return((PyArrayObject *)flux);
}

/*

   Computes the complete elliptical integral of the third kind using
   the algorithm of Bulirsch (1965):

   Bulirsch 1965, Numerische Mathematik, 7, 78
   Bulirsch 1965, Numerische Mathematik, 7, 353

 INPUTS:

    n,k - int(dtheta/((1-n*sin(theta)^2)*sqrt(1-k^2*sin(theta)^2)),0, pi/2)

 RESULT:

    The complete elliptical integral of the third kind

 -- translated from the ellpic_bulirsch.pro routine from EXOFAST
    (Eastman et al. 2013, PASP 125, 83) by Laura Kreidberg (7/22/15)
*/

inline double ellpic_bulirsch(double n, double k)
{
    double kc = sqrt(1.-k*k);
    double p = sqrt(n + 1.);
    double m0 = 1.;
    double c = 1.;
    double d = 1./p;
    double e = kc;
    double f, g;

    int nit = 0;

    while(nit < 10000)
    {
        f = c;
        c = d/p + c;
        g = e/p;
        d = 2.*(f*g + d);
        p = g + p;
        g = m0;
        m0 = kc + m0;
        if(fabs(1.-kc/g) > 1.0e-8)
        {
            kc = 2.*sqrt(e);
            e = kc*m0;
        }
        else
        {
            return 0.5*M_PI*(c*m0+d)/(m0*(m0+p));
        }
        nit++;
    }

    #ifndef _OPENACC
    printf("Convergence failure in ellpic_bulirsch\n");
    #endif
    return 0;
}

inline double ellec(double k)
{
    double m1, a1, a2, a3, a4, b1, b2, b3, b4, ee1, ee2, ellec;
    // Computes polynomial approximation for the complete elliptic
    // integral of the first kind (Hasting's approximation):
    m1 = 1.0 - k*k;
    a1 = 0.44325141463;
    a2 = 0.06260601220;
    a3 = 0.04757383546;
    a4 = 0.01736506451;
    b1 = 0.24998368310;
    b2 = 0.09200180037;
    b3 = 0.04069697526;
    b4 = 0.00526449639;
    ee1 = 1.0 + m1*(a1 + m1*(a2 + m1*(a3 + m1*a4)));
    ee2 = m1*(b1 + m1*(b2 + m1*(b3 + m1*b4)))*log(1.0/m1);
    ellec = ee1 + ee2;
    return ellec;
}

inline double ellk(double k)
{
    double a0, a1, a2, a3, a4, b0, b1, b2, b3, b4, ellk,  ek1, ek2, m1;
    // Computes polynomial approximation for the complete elliptic
    // integral of the second kind (Hasting's approximation):
    m1 = 1.0 - k*k;
    a0 = 1.38629436112;
    a1 = 0.09666344259;
    a2 = 0.03590092383;
    a3 = 0.03742563713;
    a4 = 0.01451196212;
    b0 = 0.5;
    b1 = 0.12498593597;
    b2 = 0.06880248576;
    b3 = 0.03328355346;
    b4 = 0.00441787012;
    ek1 = a0 + m1*(a1 + m1*(a2 + m1*(a3 + m1*a4)));
    ek2 = (b0 + m1*(b1 + m1*(b2 + m1*(b3 + m1*b4))))*log(m1);
    ellk = ek1 - ek2;
    return ellk;
}

static char _quadratic_ld_doc[] = "This extension module returns a limb darkened light curve for a quadratic stellar intensity profile.";

static PyMethodDef _quadratic_ld_methods[] = {
  {"_quadratic_ld", _quadratic_ld, METH_VARARGS, _quadratic_ld_doc},{NULL}};

#if PY_MAJOR_VERSION >= 3
    static struct PyModuleDef _quadratic_ld_module = {
        PyModuleDef_HEAD_INIT,
        "_quadratic_ld",
        _quadratic_ld_doc,
        -1,
        _quadratic_ld_methods
    };

    PyMODINIT_FUNC
    PyInit__quadratic_ld(void)
    {
        PyObject* module = PyModule_Create(&_quadratic_ld_module);
        if(!module)
        {
            return NULL;
        }
        import_array();
        return module;
    }
#else

    void init_quadratic_ld(void)
    {
          Py_InitModule("_quadratic_ld", _quadratic_ld_methods);
        import_array();
    }
#endif

