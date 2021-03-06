"""
Here we fit to the darksy data.
"""
import aemHMF
import numpy as np
import matplotlib.pyplot as plt
plt.rc("text", usetex=True)
plt.rc("font", size=18, family="serif")

Volume = 8000.**3 #(Mpc/h)^3

#Scale factors that darksky is at
scale_factors = np.array([0.3333, 0.5, 0.6667, 0.8, 1.0])

def get_darksky_data(a): #Masses are Msun/h, volume is [Mpc^3/h^3] comoving
    #path = "darksky_hmfs/ds14_a_halos_%.4f.hist8_m200b"%a
    #bin_center_mass, dndlM, sigma, dlogsdlogm, lower_pmass, n, expected, dm, ds, dlnm, dlns = np.genfromtxt(path, skip_header=12, unpack=True)
    lMlo, lMhi, N, Mave = np.loadtxt("othersims_hmfs/darksky_MF_a1.0.txt", unpack=True)
    C = np.loadtxt("othersims_hmfs/darksky_C_a1.0.txt")
    return lMlo, lMhi, N, Mave, C

def get_cc_prediction(bin_center_mass, a, cosmo):
    import cosmocalc as cc
    cc.set_cosmology(cosmo)
    return np.array([cc.tinker2008_mass_function(M, a, 200)*M for M in bin_center_mass])

if __name__ == "__main__":
    a = 1.0 #this is the scale factor of this darksky file
    z = 1./a - 1.
    lMlo, lMhi, N, Mave, C = get_darksky_data(a)
    Mlo = 10**lMlo
    Mhi = 10** lMhi
    M_bins = np.array([Mlo, Mhi]).T
    edges = np.append(Mlo, Mhi[-1])
    err = np.sqrt(np.diag(C))

    Ol = 0.7048737821671822
    Om = 0.295037918703847
    h = 0.6880620000000001
    Ob = 0.04676431995034128
    sig8 = 0.8344
    w0 = -1.0
    ns = 0.9676
    cosmo = {"om":Om, "ob":Ob, "ol":Ol, "h":h, "s8":sig8, "ns":ns, "w0":w0, "wa":0.0, "Neff":3.046}
    hmf = aemHMF.Aemulus_HMF()
    hmf.set_cosmology(cosmo)
    N_aem = hmf.n_bins(M_bins, a)*Volume
    Nreal = 30
    fs = hmf.residual_realization(Mave, a, Nreal)

    fig, axarr = plt.subplots(2, sharex=True)
    axarr[0].errorbar(Mave, N, err, ls='', marker='o', c='k', label="Darksky")
    axarr[0].loglog(Mave, N_aem, ls='-', c='b', label=r"Aemulus")
    #axarr[0].loglog(M, dndlM_cc, ls='-', c='r', label="Tinker08")
    axarr[0].set_ylabel(r"Number in bin")
    axarr[0].legend(loc=0, frameon=False)

    pdiff = (N - N_aem)/N_aem
    pde = err/N_aem
    
    axarr[1].errorbar(Mave, pdiff, pde, c='b', ls='-')
    #axarr[1].plot(M, pdcc, c='r', ls='-')
    axarr[1].set_ylabel(r"$\frac{N-N_{aem}}{N_{aem}}$")
    axarr[1].set_xlabel(r"$M [{\rm M}_\odot h^{-1}]$")
    ylim = .1
    axarr[1].set_ylim(-ylim, ylim)
    axarr[0].set_title(r"z=%.2f"%z)
    axarr[1].axhline(0, c='k', ls='--')

    for i in range(Nreal):
        Ns = N_aem * (1+fs[i])
        #axarr[0].loglog(Mave, Ns, ls='-', c='r', alpha=0.2)
        axarr[1].plot(Mave, (N-Ns)/N_aem, c="b", alpha=0.2)
    
    plt.subplots_adjust(hspace=0, left=0.15, bottom=0.15)
    plt.show()
