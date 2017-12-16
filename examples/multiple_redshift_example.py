import aemHMF
import numpy as np
import matplotlib.pyplot as plt
import aemulus_data as AD
plt.rc("text", usetex=True)

#Rotated truth
rtruth = np.array([-2.339125542802287516e-01,1.677880868088988064e+00,-2.970725392927930919e-01,-1.811569826499570546e+00,-8.172432603742930679e-01,4.115225664991982168e-01,5.014313130135825691e-01,6.339511118396800038e-01])

#Truths seem to not be working...
truth = np.array([2.129564190209325858e+00,2.723050117494064670e-01,1.030583349035203877e+00,2.615788704277006715e-01,4.399726014100137084e-01,1.116731634545981983e-01,1.302959111993434504e+00,2.951828692688652223e-01])


if __name__ == "__main__":
    Volume = 1050.**3 #Mpc^3/h^3
    box = 0
    Obh2, Och2, w, ns, ln10As, H0, Neff, sig8 = AD.get_building_box_cosmologies()[box]
    cosmo = {"Obh2":Obh2, "Och2":Och2, "H0":H0, "ln10^{10}A_s":ln10As, "n_s":ns, "w0":w, "N_eff":Neff}

    hmf = aemHMF.Aemulus_HMF()
    hmf.set_cosmology(cosmo)
    print hmf.tinkerMF.t08_slopes_intercepts
    hmf.tinkerMF.t08_slopes_intercepts = truth
    
    f, axarr = plt.subplots(2, sharex=True)
    sfs = AD.get_scale_factors()
    zs = 1./sfs - 1
    colors = [plt.get_cmap("seismic")(ci) for ci in np.linspace(1.0, 0.0, len(sfs))]
    for snapshot in range(len(sfs)):
        if snapshot < 0: continue
        z = zs[snapshot]
        lMlo, lMhi, N, Mtot = AD.get_building_box_binned_mass_function(box, snapshot).T
        M_bins = 10**np.array([lMlo, lMhi]).T
        M = Mtot/N
        cov = AD.get_building_box_binned_mass_function_covariance(box, snapshot)
        err = np.sqrt(np.diag(cov))

        N_aem = hmf.n_bins(M_bins, z)*Volume
        pdiff = (N-N_aem)/N_aem
        pdiff_err = err/N_aem

        icov = np.linalg.inv(cov)
        X = N - N_aem
        chi2 = np.dot(X, np.dot(icov, X))
        print "chi2(%d) = %.2f"%(snapshot, chi2)
    
        axarr[0].errorbar(M, N, err, ls='', marker='.', c=colors[snapshot], label=r"$z=%.2f$"%z)
        axarr[0].loglog(M, N_aem, ls='-', c=colors[snapshot])
        axarr[1].errorbar(M, pdiff, pdiff_err, c=colors[snapshot], ls='-', marker='.', markersize=1)
    axarr[1].axhline(0, c='k', ls='--')
    axarr[0].set_yscale('log')
    axarr[0].set_ylim(1, 1e6)
    plim = 0.07
    axarr[1].set_ylim(-plim, plim)
    axarr[1].set_ylim(-.1, 0.7)
    plt.subplots_adjust(hspace=0, wspace=0)
    axarr[0].legend(loc=0, frameon=0, fontsize=8)
    axarr[1].set_xlabel(r"${\rm Mass}\ [{\rm M_\odot}/h]$")
    axarr[1].set_ylabel(r"\% Diff")
    axarr[0].set_ylabel(r"$N$")
    plt.show()
