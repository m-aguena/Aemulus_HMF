from cosmosis.datablock import names, option_section
from scipy.interpolate import RectBivariateSpline

import sys
sys.path.append('/home/aguena/git_codes/Aemulus_HMF/aemHMF')
from tinkerMF import tinkerMF

import numpy as np

def setup(options):

    # We use a class as a namespace here to return config information

    class config:

        # Name of the section to get matter power from
        sigma_dir  = options.get_string(option_section, "sigma")

    return config

def execute(block, config):
    '''
    CAMB parameters:
        omega_b:  real, baryon density fraction today
        omega_lambda:  real, dark energy density fraction today
        omega_k:  real, curvature density fraction today (default 0.0)
        hubble:  real, hubble parameter H0 (km/s/Mpc)
        A_s:  real, scalar spectrum primordial amplitude (ignored in background/thermal mode)
        n_s:  real, scalar spectral index (ignored in background/thermal mode)
        w:  real, w(z=0) equation of state of dark energy (default -1.0). ignored if use_tabulated_w=T
        wa:  real, equation of state parameter w(z) = w_0 + w_a z / (1+z)  (default 0.0). ignored if use_tabulated_w=T
        massless_nu:  real, effective number of massless neutrinos (default 3.046)
    '''

    cosm = {"om":0.3, "ob":0.05, "ol":1.-0.3, "ok":0.0, "h":0.7,
             "s8":0.77, "ns":0.96, "w0":-1.0, "wa":0.0, "Neff":3.0}

    cosm_dic = {
        # used by emulator
        "om": "omega_m",
        "ob": "omega_b",
        "w0": "w",
        "ns": "n_s",
        "h": "hubble",
        #"Neff": "massless_nu",
        #"s8": "A_s",# still has to be completed
        # not used by emulator
        "ol": "omega_lambda",
        "ok": "omega_k",
        "wa": "wa",
            }
    for n, p in cosm_dic.items():
        cosm[n] = block[names.cosmological_parameters, p]
    cosm['h'] *= 0.01

    print('cosm:',cosm)
        
    # Just a simple rename for clarity.
    sigma_dir = config.sigma_dir

    # Load sigma(M) from the block
    #m, z, sig2 = block.get_grid(sigma_dir, "m", "z", "sigma2")
    #z, m, sig2 = block.get_grid(sigma_dir, "z", "m", "sigma2")
    m = block[sigma_dir, 'm']
    z = block[sigma_dir, 'z']
    sig = np.sqrt(block[sigma_dir, 'sigma2'])

    sig_interp = RectBivariateSpline(m, z, sig)

    class CosmosisSigma():
        def set_cosmology(self, arg):
            pass
        def sigmaMtophat(self, M, a):
            z = 1./a - 1.
            return sig_interp(np.log10(M), z)[0][0]

    '''
    a = 1./(1.+z)
    aem = tinkerMF()
    aem.set_cosmology(cosm)
    mf = np.array([[aem.dndlM(m_, a_) for a_ in a] for m_ in 10**m])
    '''

    aem = tinkerMF(sigma_obj=CosmosisSigma())
    aem.set_cosmology(cosm)
    mf = np.array([[aem.dndlM(m_, a_) for a_ in 1./(1.+z)] for m_ in 10**m])

    #block.put_grid("mf", "mf", mf, "z", z_vec, "sigma2", sigma_m.T)
    print mf.shape
    block["mf", "mf" ] = mf

    #We tell CosmoSIS that everything went fine by returning zero
    return 0

def cleanup(config):
    # Usually python modules do not need to do anything here.
    # We just leave it in out of pedantic completeness.
    pass


