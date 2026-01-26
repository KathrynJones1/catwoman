# catwoman: a batman extension to generate morning/evening terminator transit lightcurves
# Copyright (C) 2019 Kathryn Jones & Néstor Espinoza
#
# This program incorporates a modified version of the batman package: fast computation of exoplanet transit light curves
# Copyright (C) 2015 Laura Kreidberg
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
from . import _nonlinear_ld
from . import _quadratic_ld
from . import _logarithmic_ld
from . import _exponential_ld
from . import _power2_ld
from . import _custom_ld
from . import _rsky
from . import _eclipse
from math import pi
import multiprocessing
from . import openmp

__all__ = ['TransitModel', 'TransitParams']

def wrapper(func, *args, **kwargs):
    def wrapped():
        return func(*args, **kwargs)
    return wrapped

class TransitModel(object):
	"""
	Class for generating model asymmetric transit light curves.	

	:param params: A :attr:`TransitParams` object containing the physical parameters of the transit
	:type params: a `TransitParams` instance

	:param t: Array of times at which to calculate the model.
	:type t: ndarray 

	:param max_err: Error tolerance/truncation error (in parts per million) for the model.
	:type max_err: float, optional

	:param fac: Scale factor for integration step size
	:type fac: float, optional

	:param supersample_factor: Number of points subdividing exposure
	:type supersample_factor: integer, optional

	:param exp_time: Exposure time (in same units as `t`)
	:type exp_time: double, optional

	:Example:
	
	>>> model = catwoman.TransitModel(params, t, max_err = 1.0)
	"""
	#:param nthreads: Number of threads to use for parallelization.
        #:type nthreads: int, optional

        #:param transittype: Type of transit ("primary" or "secondary")
        #:type transittype: string, optional

	def __init__(self, params, t, max_err=1.0, nthreads = 1, fac = None, transittype = "primary", supersample_factor = 1, exp_time = 0.):
		#checking for invalid input
		if  (params.limb_dark == "uniform" and params.u!=None) or (params.limb_dark == "linear" and len(params.u) != 1) or \
		    (params.limb_dark == "quadratic" and len(params.u) != 2) or (params.limb_dark == "logarithmic" and len(params.u) != 2) or \
		    (params.limb_dark == "exponential" and len(params.u) != 2) or (params.limb_dark == "squareroot" and len(params.u) != 2) or \
		    (params.limb_dark == "power2" and len(params.u) != 2) or \
		    (params.limb_dark == "nonlinear" and len(params.u) != 4):
			raise Exception("Incorrect number of coefficients for " +params.limb_dark + " limb darkening; u should have the form:\n \
			 u = [] for uniform LD\n \
			 u = [u1] for linear LD\n \
  			 u = [u1, u2] for quadratic, logarithmic, exponential, squareroot, and power2 LD\n \
			 u = [u1, u2, u3, u4] for nonlinear LD, or\n \
		         u = [u1, ..., un] for custom LD") 
		if params.limb_dark not in ["uniform", "linear", "quadratic", "logarithmic", "exponential", "squareroot", "nonlinear", "power2", "custom"]: 
			raise Exception("\""+params.limb_dark+"\""+" limb darkening not supported; allowed options are:\n \
				uniform, linear, quadratic, logarithmic, exponential, squareroot, nonlinear, power2, custom")
		if max_err < 0.001: raise Exception("The lowest allowed value for max_err is 0.001. For more accurate calculation, set the integration step size explicitly with the fac parameter.")
		if transittype not in ["primary", "secondary"]: raise Exception("Allowed transit types are \"primary\" and \"secondary\".")
		if (supersample_factor > 1 and exp_time <= 0.): raise Exception("Please enter a valid exposure time (exp_time must be greater than 0 to calculate super-sampled light curves).")
		if (not isinstance(t, np.ndarray)): raise Exception("Times t must be a numpy array (not a list).")
		if (params.phi>90) or (params.phi<(-90)): raise Exception("Please enter a valid phi angle (params.phi must be between 90 and -90 degrees).")  
		
		#initializes model parameters
		self.t = t
		self.t0 = params.t0
		self.per = params.per
		self.rp = params.rp
		self.rp2 = params.rp2 #added this
		self.a = params.a
		self.inc = params.inc
		self.ecc = params.ecc
		self.w = params.w
		self.u = params.u
		self.limb_dark = params.limb_dark
		self.fp = params.fp
		self.t_secondary = params.t_secondary
		self.max_err = max_err
		self.supersample_factor = supersample_factor
		self.exp_time = exp_time
		self.inverse = False
		self.twocircles = False        
		self.phi = params.phi*pi/180          #convert phi from degrees to radians
		self.b = params.a*np.cos(params.inc*pi/180)*((1-params.ecc*params.ecc)/(1-params.ecc*np.sin(params.w*pi/180)))            	      		

		if self.supersample_factor > 1:  # IJMC: now do it quicker, with no loops:
			t_offsets = np.linspace(-self.exp_time/2., self.exp_time/2., self.supersample_factor)
			self.t_supersample = (t_offsets + self.t.reshape(self.t.size, 1)).flatten()
			self.t = self.t_supersample
		else: self.t_supersample = self.t
		
		self.phi2=np.zeros(len(self.t))
		self.Y=np.zeros(len(self.t))
		self.X=np.zeros(len(self.t))
		self.psi=np.zeros(len(self.t))				
		self.big_vector = np.zeros(4*len(self.t))
	
		#Checking if there are two semi-circles	
		if self.rp2==None:
			self.twocircles=False
			self.rp2=0.0
		else:
			self.twocircles=True


		#handles the case of inverse transits (rp < 0)
		if self.rp < 0 : 
			self.rp = -1.*self.rp
			params.rp = -1.*params.rp
			self.inverse = True

		if transittype == "primary": self.transittype = 1
		else: 
			self.transittype = 2
			params.t0 = self.get_t_conjunction(params)		
		
		#Setting a fac if uniform limb darkening or if already specified, if not calculate it
		if fac != None: self.fac = fac
		elif self.limb_dark != "uniform": self.fac = self._get_fac()
		else: self.fac = 0.5
		
		if nthreads==None or nthreads == 1: self.nthreads=1
		else:
			if nthreads <= multiprocessing.cpu_count()and nthreads >1 and openmp.detect(): self.nthreads = nthreads
			else: 
				if nthreads > multiprocessing.cpu_count(): raise Exception("Maximum number of threads is "+'{0:d}'.format(multiprocessing.cpu_count()))
				elif nthreads <= 1: raise Exception("Number of threads must be between 2 and {0:d}".format(multiprocessing.cpu_count()))
				else: raise Exception("OpenMP not enabled: do not set the nthreads parameter")
		self.big_vector = _rsky._rsky(self.t_supersample, params.t0, params.per, params.a, params.inc*pi/180., params.ecc, params.w*pi/180., self.transittype, self.nthreads)
		
		self.ds=self.big_vector[0:int(len(self.big_vector)/4)]
		self.Y=self.big_vector[int(len(self.big_vector)/4):int(len(self.big_vector)*2/4)]
		self.psi=self.big_vector[int(len(self.big_vector)*2/4):int(len(self.big_vector)*3/4)]
		self.X=self.big_vector[int(len(self.big_vector)*3/4):int(len(self.big_vector))]		
		
		#Correcting phi angle to orbital motion
		self.phi2 = self.phi + self.psi

	def calc_err(self, plot = False):
		"""

		Calculate maximum error for a specific transit light curve calculation.
			
		:param plot: If ``True``, plots the error in the light curve model against the separation of centres, d.
		:type plot: bool

		:return: Truncation error (in parts per million)
		:rtype: float

		"""
		if (self.limb_dark in ["logarithmic", "exponential", "nonlinear", "squareroot", "power2", "custom", "quadratic", "linear", "uniform"]):
			
			ds = np.linspace(0., 1.1, 500)
			fac_lo = 5.0e-4
			bfac = np.zeros(500)
			phifac=np.asarray([self.phi for i in range(0,500)])

			if self.limb_dark == "nonlinear":
				f0 = _nonlinear_ld._nonlinear_ld(ds, self.rp, self.u[0], self.u[1], self.u[2], self.u[3], fac_lo, self.nthreads, phifac, bfac, 0, self.rp2, self.twocircles)
				f = _nonlinear_ld._nonlinear_ld(ds, self.rp, self.u[0], self.u[1], self.u[2], self.u[3], self.fac, self.nthreads, phifac, bfac, 0, self.rp2, self.twocircles)
			elif self.limb_dark == "squareroot":
				f0 = _nonlinear_ld._nonlinear_ld(ds, self.rp, self.u[1], self.u[0], 0., 0., fac_lo, self.nthreads, phifac, bfac, 0, self.rp2,self.twocircles)
				f = _nonlinear_ld._nonlinear_ld(ds, self.rp, self.u[1], self.u[0], 0., 0., self.fac, self.nthreads, phifac, bfac, 0, self.rp2,self.twocircles)
			elif self.limb_dark == "exponential":
				f0 = _exponential_ld._exponential_ld(ds, self.rp, self.u[0], self.u[1], fac_lo, self.nthreads, phifac, bfac, 0, self.rp2,self.twocircles)
				f = _exponential_ld._exponential_ld(ds, self.rp, self.u[0], self.u[1], self.fac, self.nthreads, phifac, bfac, 0, self.rp2,self.twocircles)
			elif self.limb_dark == "logarithmic":
				f0 = _logarithmic_ld._logarithmic_ld(ds, self.rp, self.u[0], self.u[1], fac_lo, self.nthreads, phifac, bfac, 0, self.rp2,self.twocircles)
				f = _logarithmic_ld._logarithmic_ld(ds, self.rp, self.u[0], self.u[1], self.fac, self.nthreads, phifac, bfac, 0, self.rp2,self.twocircles)
			elif self.limb_dark == "power2":
				f0 = _power2_ld._power2_ld(ds, self.rp, self.u[0], self.u[1], fac_lo, self.nthreads, phifac, bfac, 0, self.rp2,self.twocircles)
				f = _power2_ld._power2_ld(ds, self.rp, self.u[0], self.u[1], self.fac, self.nthreads, phifac, bfac, 0, self.rp2,self.twocircles)
			elif self.limb_dark == "quadratic":
                                f0 = _quadratic_ld._quadratic_ld(ds, self.rp, self.u[0], self.u[1], fac_lo, self.nthreads, phifac, bfac, 0, self.rp2, self.twocircles)
                                f = _quadratic_ld._quadratic_ld(ds, self.rp, self.u[0], self.u[1], self.fac, self.nthreads, phifac, bfac, 0, self.rp2, self.twocircles)
			elif self.limb_dark == "linear":
                                f0 = _quadratic_ld._quadratic_ld(ds, self.rp, self.u[0], 0., fac_lo, self.nthreads, phifac, bfac, 0, self.rp2, self.twocircles)
                                f = _quadratic_ld._quadratic_ld(ds, self.rp, self.u[0], 0., self.fac, self.nthreads, phifac, bfac, 0, self.rp2, self.twocircles)
			elif self.limb_dark == "uniform":
                                f0 = _quadratic_ld._quadratic_ld(ds, self.rp, 0., 0., fac_lo, self.nthreads, phifac, bfac, 0, self.rp2, self.twocircles)
                                f = _quadratic_ld._quadratic_ld(ds, self.rp, 0., 0., self.fac, self.nthreads, phifac, bfac, 0, self.rp2, self.twocircles)
			else:
				f0 = _custom_ld._custom_ld(ds, self.rp, self.u[0], self.u[1], self.u[2], self.u[3], self.u[4], self.u[5], fac_lo, self.nthreads, phifac, bfac, 0, self.rp2,self.twocircles)
				f =  _custom_ld._custom_ld(ds, self.rp, self.u[0], self.u[1], self.u[2], self.u[3], self.u[4], self.u[5], self.fac, self.nthreads, phifac, bfac, 0, self.rp2,self.twocircles)
	
			err = np.max(np.abs(f-f0))*1.0e6
			if plot == True:
				import matplotlib.pyplot as plt
				plt.plot(ds, 1.0e6*(f-f0), color='k')
				plt.xlabel("d (separation of centers)")
				plt.ylabel("Error (ppm)") 
				plt.show()

			return err
		else: raise Exception("Function calc_err not valid for " + self.limb_dark + " limb darkening")

	def _get_fac(self):
		if (self.limb_dark in ["logarithmic", "exponential", "nonlinear", "squareroot", "power2", "custom", "quadratic", "linear", "uniform"]):
			nthreads = 1
			fac_lo, fac_hi = 5.0e-4, 1.
			ds = np.linspace(0., np.maximum(1.+self.rp,1.+self.rp2), 1000)
			bfac=np.zeros(1000)
			phifac=np.asarray([self.phi for i in range(0,1000)])

			if self.limb_dark == "nonlinear": f0 = _nonlinear_ld._nonlinear_ld(ds, self.rp, self.u[0], self.u[1], self.u[2], self.u[3], fac_lo, nthreads, phifac, bfac, ds, self.rp2,self.twocircles)
			elif self.limb_dark == "squareroot": f0 = _nonlinear_ld._nonlinear_ld(ds, self.rp, self.u[1], self.u[0], 0., 0., fac_lo, nthreads, phifac, bfac, ds, self.rp2,self.twocircles)
			elif self.limb_dark == "exponential": f0 = _exponential_ld._exponential_ld(ds, self.rp, self.u[0], self.u[1], fac_lo, nthreads, phifac, bfac, ds, self.rp2,self.twocircles)
			elif self.limb_dark == "logarithmic": f0 = _logarithmic_ld._logarithmic_ld(ds, self.rp, self.u[0], self.u[1], fac_lo, nthreads, phifac, bfac, ds, self.rp2,self.twocircles)
			elif self.limb_dark == "power2": f0 = _power2_ld._power2_ld(ds, self.rp, self.u[0], self.u[1], fac_lo, nthreads, phifac, bfac, ds, self.rp2, self.twocircles)
			elif self.limb_dark == "quadratic": f0 = _quadratic_ld._quadratic_ld(ds, self.rp, self.u[0], self.u[1], fac_lo, nthreads, phifac, bfac, ds, self.rp2, self.twocircles)
			elif self.limb_dark == "linear": f0 = _quadratic_ld._quadratic_ld(ds, self.rp, self.u[0], 0., fac_lo, nthreads, phifac, bfac, ds, self.rp2, self.twocircles)
			elif self.limb_dark == "uniform": f0 = _quadratic_ld._quadratic_ld(ds, self.rp, 0., 0., fac_lo, nthreads, phifac, bfac, ds, self.rp2, self.twocircles)
			else: f0 = _custom_ld._custom_ld(ds, self.rp, self.u[0], self.u[1], self.u[2], self.u[3], self.u[4], self.u[5], fac_lo, nthreads, phifac, bfac, ds, self.rp2,self.twocircles)

			n = 0
			err = 0.
			while(err > self.max_err or err < 0.99*self.max_err):
				#fac = (fac_lo + fac_hi)/2. #arithmetic mean
				fac = np.sqrt(fac_lo*fac_hi) #geometric mean(arithmetic mean in log space)

				if self.limb_dark == "nonlinear": f = _nonlinear_ld._nonlinear_ld(ds, self.rp, self.u[0], self.u[1], self.u[2], self.u[3], fac, nthreads, phifac, bfac, ds ,self.rp2,self.twocircles)
				elif self.limb_dark == "squareroot": f = _nonlinear_ld._nonlinear_ld(ds, self.rp, self.u[1], self.u[0], 0., 0., fac, nthreads, phifac, bfac, ds,self.rp2,self.twocircles)
				elif self.limb_dark == "exponential": f = _exponential_ld._exponential_ld(ds, self.rp, self.u[0], self.u[1], fac, nthreads, phifac, bfac, ds,self.rp2,self.twocircles)
				elif self.limb_dark == "logarithmic": f = _logarithmic_ld._logarithmic_ld(ds, self.rp, self.u[0], self.u[1], fac, nthreads, phifac, bfac, ds,self.rp2,self.twocircles)
				elif self.limb_dark == "power2": f = _power2_ld._power2_ld(ds, self.rp, self.u[0], self.u[1], fac, nthreads, phifac, bfac, ds, self.rp2, self.twocircles)
				elif self.limb_dark == "quadratic": f = _quadratic_ld._quadratic_ld(ds, self.rp, self.u[0], self.u[1], fac, nthreads, phifac, bfac, ds, self.rp2, self.twocircles)
				elif self.limb_dark == "linear": f = _quadratic_ld._quadratic_ld(ds, self.rp, self.u[0], 0., fac, nthreads, phifac, bfac, ds, self.rp2, self.twocircles)
				elif self.limb_dark == "uniform": f = _quadratic_ld._quadratic_ld(ds, self.rp, 0., 0., fac, nthreads, phifac, bfac, ds, self.rp2, self.twocircles)
				else: f = _custom_ld._custom_ld(ds, self.rp, self.u[0], self.u[1], self.u[2], self.u[3], self.u[4], self.u[5], fac, nthreads, phifac, bfac, ds ,self.rp2,self.twocircles)

				err = np.max(np.abs(f-f0))*1.0e6
				if err > self.max_err: fac_hi = fac	
				else: fac_lo = fac
				n += 1
				if n > 1e3: raise Exception("Convergence failure in calculation of scale factor for integration step size")
			return fac
		else: return 0.5 
		#else: return None
	
	def light_curve(self, params):
		"""
		Calculates and returns a model asymmetric light curve.

		:param params: Transit parameters
		:type params: A `TransitParams` instance

		:return: Relative flux 
		:rtype: ndarray

		:Example:

		>>> flux = model.light_curve(params)
		"""
		#recalculates rsky and fac if necessary
		new_rsky=0
		if params.t0 != self.t0 or params.per != self.per or params.a != self.a or params.inc != self.inc or params.ecc != self.ecc or params.w != self.w or params.t_secondary != self.t_secondary:
			if self.transittype == 2 and params.t_secondary != self.t_secondary:
				params.t0 = self.get_t_conjunction(params)
			self.big_vector= _rsky._rsky(self.t_supersample, params.t0, params.per, params.a, params.inc*pi/180., params.ecc, params.w*pi/180., self.transittype, self.nthreads)
			
			self.ds=self.big_vector[0:int(len(self.big_vector)/4)]
			self.Y=self.big_vector[int(len(self.big_vector)/4):int(len(self.big_vector)*2/4)]
			self.psi=self.big_vector[int(len(self.big_vector)*2/4):int(len(self.big_vector)*3/4)]
			self.X=self.big_vector[int(len(self.big_vector)*3/4):int(len(self.big_vector))]
			new_rsky = 1
		if params.limb_dark != self.limb_dark: self.fac = self._get_fac()

		#updates transit params
		self.t0 = params.t0
		self.per = params.per
		self.rp = params.rp
		self.rp2 = params.rp2
		self.a = params.a
		self.inc = params.inc
		self.ecc = params.ecc
		self.w = params.w
		self.u = params.u
		self.limb_dark = params.limb_dark
		self.fp = params.fp
		self.t_secondary = params.t_secondary
		self.inverse = False
		self.twocircles = False
		self.b = params.a*np.cos(params.inc*pi/180)*((1-params.ecc*params.ecc)/(1-params.ecc*np.sin(params.w*pi/180)))
		
		#updating self.phi2 if self.phi or self.psi has changed
		if ((params.phi*pi/180)!=self.phi or (new_rsky==1)):
			self.phi = params.phi*pi/180
              
			#Correcting phi angle to orbital motion                	
			self.phi2 = self.phi + self.psi
				
		#checking again if there has been a value of rp2 supplied
		if self.rp2 == None: 
			self.twocircles = False
			self.rp2=0.0
		else: self.twocircles = True

		#handles the case of inverse transits (rp < 0)
		if self.rp < 0.: 
			self.rp = -1.*self.rp
			params.rp = -1.*params.rp
			self.inverse = True
		
		
		if self.transittype == 1:
			if params.limb_dark != self.limb_dark: raise Exception("Need to reinitialize model in order to change limb darkening option")
			if self.limb_dark == "quadratic": lc = _quadratic_ld._quadratic_ld(self.ds, params.rp, params.u[0], params.u[1], self.fac, self.nthreads, self.phi2, self.Y, self.X, self.rp2, self.twocircles)
			elif self.limb_dark == "uniform": lc = _quadratic_ld._quadratic_ld(self.ds, params.rp, 0., 0., self.fac, self.nthreads, self.phi2, self.Y, self.X, self.rp2, self.twocircles)
			elif self.limb_dark == "linear": lc = _quadratic_ld._quadratic_ld(self.ds, params.rp, params.u[0], 0.,self.fac, self.nthreads, self.phi2, self.Y, self.X, self.rp2, self.twocircles)
			elif self.limb_dark == "nonlinear": lc = _nonlinear_ld._nonlinear_ld(self.ds, params.rp, params.u[0], params.u[1], params.u[2], params.u[3], self.fac, self.nthreads, self.phi2, self.Y, self.X,self.rp2,self.twocircles)
			elif self.limb_dark == "squareroot": lc = _nonlinear_ld._nonlinear_ld(self.ds, params.rp, params.u[1], params.u[0], 0., 0., self.fac, self.nthreads, self.phi2, self.Y, self.X,self.rp2,self.twocircles)
			elif self.limb_dark == "logarithmic": lc = _logarithmic_ld._logarithmic_ld(self.ds, params.rp, params.u[0], params.u[1], self.fac, self.nthreads, self.phi2, self.Y, self.X,self.rp2,self.twocircles)
			elif self.limb_dark == "exponential": lc = _exponential_ld._exponential_ld(self.ds, params.rp, params.u[0], params.u[1], self.fac, self.nthreads, self.phi2, self.Y, self.X,self.rp2,self.twocircles)
			elif self.limb_dark == "power2": lc = _power2_ld._power2_ld(self.ds, params.rp, params.u[0], params.u[1], self.fac, self.nthreads, self.phi2, self.Y, self.X,self.rp2,self.twocircles)
			elif self.limb_dark == "custom": lc = _custom_ld._custom_ld(self.ds, params.rp, params.u[0], params.u[1], params.u[2], params.u[3], params.u[4], params.u[5], self.fac, self.nthreads, self.phi2, self.Y, self.X,self.rp2,self.twocircles)
			else: raise Exception("Invalid limb darkening option")

			if self.inverse == True: lc = 2. - lc

		else: lc = _eclipse._eclipse(self.ds, params.rp, params.fp, self.nthreads)			
		if self.supersample_factor == 1: return lc
		else: return np.mean(lc.reshape(-1, self.supersample_factor), axis=1)

	def _get_phase(self, params, position):
		if position == "periastron": TA = 0.
		elif position == "primary": TA = pi/2. - params.w*pi/180.
		elif position == "secondary": TA = 3.*pi/2. - params.w*pi/180.
		
		E = 2.*np.arctan(np.sqrt((1. - params.ecc)/(1. + params.ecc))*np.tan(TA/2.))
		M = E - params.ecc*np.sin(E)
		return M/2./pi
	
	def get_t_periastron(self, params):
		"""
		Return the time of periastron passage (calculated using `params.t0`).
		"""
		phase = self._get_phase(params, "primary")
		return params.t0 - params.per*phase

	def get_t_secondary(self, params):
		"""
		Return the time of secondary eclipse center (calculated using `params.t0`).
		"""
		phase = self._get_phase(params, "primary")
		phase2 = self._get_phase(params, "secondary")
		return params.t0 + params.per*(phase2-phase)

	def get_t_conjunction(self, params):
		"""
		Return the time of primary transit center (calculated using `params.t_secondary`).
		"""
		phase = self._get_phase(params, "primary")
		phase2 = self._get_phase(params, "secondary")
		return params.t_secondary + params.per*(phase-phase2)

	def get_true_anomaly(self):
		"""
		Return the true anomaly at each time
		"""
		big_vec_getf = _rsky._getf(self.t_supersample, self.t0, self.per, self.a,
							  self.inc*pi/180., self.ecc, self.w*pi/180.,
							  self.transittype, self.nthreads)
		self.f = big_vec_getf[0:len(big_vec_getf)/4]
		return self.f

class TransitParams(object):
	"""
	Object to store the physical parameters of the transit.

	:param t0: Time of inferior conjunction. 
	:type t0: float 

	:param t_secondary: Time of secondary eclipse center.
	:type t_secondary: float, optional 

	:param per: Orbital period.
	:type per: float

	:param rp: Planet radius 1 (in stellar radii).
	:type rp: float

	:param rp2: Planet radius 2 (in stellar radii) 
	:type rp2: float                               

	:param a: Semi-major axis (in stellar radii).
	:type a: float

	:param inc: Orbital inclination (in degrees).
	:type inc: float

	:param ecc: Orbital eccentricity.
	:type ecc: float

	:param w: Argument of periapse (in degrees)
	:type w: float

	:param u: List of limb darkening coefficients.
	:type u: array_like 

	:param limb_dark: Limb darkening model (choice of "nonlinear", "quadratic", "exponential", "logarithmic", "squareroot", "linear", "uniform", "power2", or "custom").
	:type limb_dark: str

	:param fp: Planet-to-star flux ratio (for secondary eclipse models).
	:type fp: float, optional

	:param phi: Rotation angle of planet, relative to horizontal (in degrees)
	:type phi: float

	.. note::  
		- Units for the orbital period and ephemeris can be anything as long as they are consistent (e.g. both in days). 

	:Example:
	
	>>> import catwoman
	>>> params = catwoman.TransitParams()
	>>> params.t0 = 0. 				#time of inferior conjunction
	>>> params.per = 1.				#orbital period	
	>>> params.rp = 0.1				#planet radius 1 (in units of stellar radii)
	>>> params.rp2 = 0.1005				#planet radius 2 (in units of stellar radii)
	>>> params.a = 15.				#semi-major axis (in units of stellar radii)
	>>> params.inc = 87.				#orbital inclination (in degrees)	
	>>> params.ecc = 0.				#eccentricity	
	>>> params.w = 90.				#longitude of periastron (in degrees) 
	>>> params.u = [0.1, 0.3] 	      	        #limb darkening coefficients
	>>> params.limb_dark = "quadratic"          	#limb darkening model
	>>> params.phi = 90.				#rotation angle of planet
	"""	
	# - The orbital path is calculated based on `t0` for primary transits and `t_secondary` for secondary eclipses.

	def __init__(self):
		self.t0 = None
		self.per = None
		self.rp = None  
		self.rp2 = None 
		self.a = None
		self.inc = None
		self.ecc = None
		self.w = None
		self.u = None
		self.limb_dark = None
		self.fp = None
		self.t_secondary = None
		self.phi = None
