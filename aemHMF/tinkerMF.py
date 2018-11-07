"""
This file contains the actual Tinker08 mass function. It uses emu to get the mass function parameters.
"""
#from classy import Class
#from cluster_toolkit import bias
import numpy as np
from scipy import special#, integrate
from scipy.interpolate import InterpolatedUnivariateSpline as IUS
#from scipy.interpolate import RectBivariateSpline as RBS
#from scipy.interpolate import interp2d as RBS
import emu

#Physical constants
G = 4.51715e-48 #Newton's gravitional constant in Mpc^3/s^2/Solar Mass
Mpcperkm = 3.240779289664256e-20 #Mpc/km; used to convert H0 to s^-1
rhocrit = 3.*(Mpcperkm*100)**2/(8*np.pi*G) #Msun h^2/Mpc^3

class tinkerMF(object):
    
    def __init__(self, tinker=False, sigma_obj=None):
        self.t_08 = emu.emu()
        self.tinker = tinker

        self.sigma_obj = sigma_obj
        if sigma_obj is None:
            import cosmocalc
            self.sigma_obj = cosmocalc


    def set_cosmology(self, cosmo_dict):
        self.cosmo_dict = cosmo_dict
        self.rhom = cosmo_dict['om']*rhocrit #Msun h^2/Mpc^3
        self.t08_slopes_intercepts = \
            self.t_08.predict_slopes_intercepts(
                    self.cos_from_dict(cosmo_dict))
        self.sigma_obj.set_cosmology(cosmo_dict)
        
    def cos_from_dict(self, cosmo_dict):
        cd = cosmo_dict
        om = cd['om']
        ob = cd['ob']
        w0 = cd['w0']
        ns = cd['ns']
        h  = cd['h']
        Neff = cd['Neff']
        s8 = cd['s8']
        H0 = h*100
        Obh2 = ob*h*h
        Och2 = (om-ob)*h*h
        return np.array([Obh2, Och2, w0, ns, H0, Neff, s8])

    def calc_normalization(self):
        #Calculates B, the normalization of the T08 MF.
        d,e,f,g = self.t08_params
        gamma_d2 = special.gamma(d*0.5)
        gamma_f2 = special.gamma(f*0.5)
        gnd2 = g**(-d*0.5)
        gnf2 = g**(-f*0.5)
        ed = e**d
        self.B = 2.0/(ed * gnd2 * gamma_d2 + gnf2 * gamma_f2)
        return
        
    def merge_t08_params(self, a):
        if self.tinker:
            d, e, f, g = 1.97, 1.00, 0.51, 1.228
        else:
            k = a-0.5
            d0,d1,f0,f1,g0,g1 = self.t08_slopes_intercepts
            d = d0 + k*d1
            e = np.array([1.0]) #Default Tinker08 value
            f = f0 + k*f1
            g = g0 + k*g1
        self.t08_params = np.array([d, e, f, g]).flatten()
        return

    def MtoR(self, M):
        return (M/(4./3.*np.pi*self.rhom))**(1./3.)/self.cosmo_dict['h'] #Lagrangian radius in Mpc

    def Mtosigma(self, M, a):
        return self.sigma_obj.sigmaMtophat(M, a)
        
    def Gsigma(self, sigma, a):
        if not hasattr(self, "a"):
            self.a = a
            self.merge_t08_params(a)
            self.calc_normalization() #Recalculate B
        if self.a != a:
            self.a = a
            self.merge_t08_params(a)
            self.calc_normalization() #Recalculate B
        d,e,f,g = self.t08_params
        return self.B*((sigma/e)**-d + sigma**-f) * np.exp(-g/sigma**2)

    def GM(self, M, a):
        if not hasattr(self, "a"):
            self.a = a
            self.merge_t08_params(a)
            self.calc_normalization() #Recalculate B
        if self.a != a:
            self.a = a
            self.merge_t08_params(a)
            self.calc_normalization() #Recalculate B
        sigma = self.Mtosigma(M, a)
        return self.Gsigma(sigma, a)

    def dndlM(self, M, a):
        gsigma = self.GM(M, a)
        dM = 1e-6*M
        dlnsiginvdm = np.log(self.Mtosigma(M-dM/2,a)/self.Mtosigma(M+dM/2, a))/dM
        return gsigma * self.rhom * dlnsiginvdm

    def n_bin(self, Mlow, Mhigh, a):
        M = np.logspace(np.log10(Mlow), np.log10(Mhigh), num=1000)
        lM = np.log(M)
        dndlM = np.array([self.dndlM(Mi, a) for Mi in M])
        spl = IUS(lM, dndlM)
        return spl.integral(np.log(Mlow), np.log(Mhigh))

    def n_bins(self, Mbins, a):
        return np.array([self.n_bin(Mbi[0], Mbi[1], a) for Mbi in Mbins])

    def peak_height(self, M, a):
        return 1.686/self.Mtosigma(M,a)

if __name__ == "__main__":
    cos = {"om":0.3,"ob":0.05,"ol":1.-0.3,"ok":0.0,"h":0.7,"s8":0.77,"ns":0.96,"w0":-1.0,"wa":0.0,"Neff":3.0}
    n = tinkerMF()
    M = np.logspace(12, 16, num=10)

    for i in xrange(0,1):
        a = 1./(1.+i)
        c = (1-a)/2
        n.set_cosmology(cos)
        dndlMcc = np.array([n.dndlM(Mi, a) for Mi in M])
