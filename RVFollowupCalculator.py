'''
Compute nRV for all TESS targets using the code cloned from
https://github.com/r-cloutier/Transiting_RVcalculator.git
'''
from imports import *
from compute_sigmaRV import *
from sigmaRV_activity import *
from sigmaRV_planets import *
from compute_nRV_GP import *
from create_pdf import *

global G
G = 6.67e-11

def nRV_calculator(Kdetsig,
                   input_planet_fname='user_planet.in',
                   input_star_fname='user_star.in',
                   input_spectrograph_fname='user_spectrograph.in',
                   input_sigRV_fname='user_sigRV.in',
                   output_fname='RVFollowupCalculator',
                   Ntrials=1):
    '''
    Compute the number of RV measurements required to detect an input 
    transiting planet around an input star with an input spectrograph at a 
    given significance.

    Parameters
    ----------
    `Kdetsig': scalar
        The desired RV semi-amplitude detection significance measured as 
        the semi-amplitude over its measurement uncertainty 
        (i.e. Kdetsig = K / sigmaK)
    `Ntrials' : scalar
        The number of nRV calculations to perform. Must be >= 1.

    Returns
    -------
    `nRV': float
        The number of RV measurements required to detect the planet's RV 
        semi-amplitude with a detection significance of `Kdetsig'
    `texp': float
        The exposure time in minutes required to achieve a SNR of SNRtarget 
        (see input_spectrograph_fname) per resolution element
    `tobs': float
        The total observing time in hours required to detect the planet's RV 
        semi-amplitude with a detection significance of `Kdetsig'

    '''    
    # get inputs
    P, rp, mp = _read_planet_input(input_planet_fname)
    mags, Ms, Rs, Teff, Z, vsini, Prot = _read_star_input(input_star_fname)
    band_strs, R, aperture, throughput, RVnoisefloor, centralwl_nm, SNRtarget, \
        transmission_threshold, texpmin, texpmax, toverhead = \
                            _read_spectrograph_input(input_spectrograph_fname)
    assert texpmin < texpmax
    assert mags.size == band_strs.size
    sigRV_phot, sigRV_act, sigRV_planet, sigRV_eff = _read_sigRV_input(input_sigRV_fname)

    # get RV noise sources if effective RV rms is not specified
    Ntrial = int(Ntrials)
    assert Ntrials > 0
    if sigRV_eff <= 0:
        logg = float(unp.nominal_values(_compute_logg(Ms, Rs)))

	if sigRV_phot > 0:
	    texp = exposure_time_calculator_per_band(mags, band_strs, aperture,
                                                     throughput, R, SNRtarget,
                                                     texpmin, texpmax)
	else:
            sigRV_phot, texp = _compute_sigRV_phot(band_strs, mags, Teff, logg, Z,
                                               	   vsini, R, aperture, throughput,
                                               	   RVnoisefloor, centralwl_nm,
                                               	   SNRtarget,
                                               	   transmission_threshold, texpmin,
                                               	   texpmax)

        sigRV_acts, sigRV_planets = np.zeros(Ntrials), np.zeros(Ntrials)
        for i in range(Ntrials):      
            sigRV_acts[i] = _get_sigRV_act() if sigRV_act < 0 else float(sigRV_act)
            sigRV_planets[i] = _get_sigRV_planets() if sigRV_planet < 0 \
                               else float(sigRV_planet)
        sigRV_effs = np.sqrt(sigRV_phot**2 + sigRV_acts**2 + sigRV_planets**2)

    else:
	sigRV_acts, sigRV_planets = np.repeat(sigRV_eff, Ntrials), np.zeros(Ntrials)
        sigRV_effs = np.repeat(sigRV_eff, Ntrials)
        texp = exposure_time_calculator_per_band(mags, band_strs, aperture,
                                                 throughput, R, SNRtarget,
                                                 texpmin, texpmax)

    # get target K measurement uncertainty
    mp = float(_get_planet_mass(rp)) if mp == 0 else float(mp)
    sigK_target = _get_sigK(Kdetsig, P, Ms, mp)

    # compute number of RVs required for a white and red noise model
    nRVs, nRVGPs = np.zeros(Ntrials), np.zeros(Ntrials)
    for i in range(Ntrials):
        GPtheta = sigRV_acts[i], Prot*3, 2., Prot, sigRV_planets[i]
        keptheta = P, rvs.RV_K(P, Ms, mp)
        nRVGPs[i] = compute_nRV_GP(GPtheta, keptheta, sigRV_phot, sigK_target)
        nRVs[i] = 2. * (sigRV_effs[i] / sigK_target)**2

    # compute total observing time
    tobss = nRVs * (texp + toverhead) / 60.  # in hours
    tobsGPs = nRVGPs * (texp + toverhead) / 60.  # in hours
    
    # write results to file
    output = [P, rp, mp, rvs.RV_K(P, Ms, mp),
              mags, Ms, Rs, Teff, Z, vsini, Prot,
              band_strs, R, aperture, throughput, RVnoisefloor, centralwl_nm*1e-3, SNRtarget,
              transmission_threshold, texpmin, texpmax, toverhead,
              sigRV_phot, sigRV_acts, sigRV_planets, sigRV_effs, sigK_target,
              texp, nRVs, nRVGPs, tobss, tobsGPs]
    #_write_results2file(output_fname, output)
    create_pdf(output_fname, output)

    return nRV, texp, tobs


def _read_planet_input(input_planet_fname):
    '''
    Read-in planetary data from the input file.
    '''
    f = open('InputFiles/%s'%input_planet_fname, 'r')
    g = f.readlines()
    f.close()
    return float(g[3]), float(g[5]), float(g[7])


def _read_star_input(input_star_fname):
    '''
    Read-in stellar data from the input file.
    '''
    f = open('InputFiles/%s'%input_star_fname, 'r')
    g = f.readlines()
    f.close()
    return np.ascontiguousarray(g[3].split(',')).astype(float), \
        float(g[5]), float(g[7]), float(g[9]), float(g[11]), float(g[13]), \
        float(g[15])


def _read_spectrograph_input(input_spectrograph_fname):
    '''
    Read-in spectrograph data from the input file.
    '''
    f = open('InputFiles/%s'%input_spectrograph_fname, 'r')
    g = f.readlines()
    f.close()
    return np.ascontiguousarray(list(g[3])[:-1]), float(g[5]), float(g[7]), \
        float(g[9]), float(g[11]), float(g[13]), float(g[15]), float(g[17]), \
        float(g[19]), float(g[21]), float(g[23])


def _read_sigRV_input(input_sigRV_fname):
    '''
    Read-in RV noise source data from the input file.
    '''
    f = open('InputFiles/%s'%input_sigRV_fname, 'r')
    g = f.readlines()
    f.close()
    return float(g[3]), float(g[5]), float(g[7]), float(g[9])


def _compute_logg(Ms, Rs):
    '''
    Compute stellar logg in cgs units.
    '''
    Ms, Rs = rvs.Msun2kg(Ms), rvs.Rsun2m(Rs)
    return unp.log10(G * Ms / Rs**2 * 1e2)


def _compute_sigRV_phot(band_strs, mags, Teff, logg, Z, vsini, R, aperture,
                        throughput, RVnoisefloor, centralwl_nm, SNRtarget,
                        transmission_threshold, texpmin, texpmax):
    '''
    Calculate the photon-noise limited RV precision over the spectrograph's 
    full spectral domain.
    '''
    # get round values for PHOENIX stellar models
    Teffs = np.append(np.arange(23e2,7e3,1e2), np.arange(7e3,121e2,2e2))
    Teff_round = Teffs[abs(Teffs-Teff) == np.min(abs(Teffs-Teff))][0]
    loggs = np.arange(0, 6.1, .5)
    logg_round = loggs[abs(loggs-logg) == np.min(abs(loggs-logg))][0]
    Zs = np.append(np.arange(-4,-1,dtype=float), np.arange(-1.5,1.5,.5))
    Z_round = Zs[abs(Zs-Z) == np.min(abs(Zs-Z))][0]

    # get exposure time
    texp = exposure_time_calculator_per_band(mags, band_strs, aperture,
                                             throughput, R, SNRtarget,
                                             texpmin=texpmin, texpmax=texpmax)

    # compute sigmaRV in each band for a fixed texp
    sigmaRVs = np.zeros(mags.size)
    for i in range(sigmaRVs.size):
        t0 = time.time()
        wl, spec = get_reduced_spectrum(Teff_round, logg_round, Z_round, vsini,
                                        band_strs[i], R, centralwl_nm*1e-3,
                                        SNRtarget)
        sigmaRVs[i] = compute_sigmaRV(wl, spec, mags[i], band_strs[i], texp,
                                      aperture, throughput, R,
                                      transmission_threshold, SNRtarget)
        print 'Took %.1f seconds\n'%(time.time()-t0)
        
    # compute sigmaRV over all bands
    sigRV_phot = 1 / np.sqrt(np.sum(1. / sigmaRVs**2))
    sigRV_phot = sigRV_phot if sigRV_phot > RVnoisefloor \
                 else float(RV_noisefloor)
    return sigRV_phot, texp


def _get_planet_mass(rps, Fs=336.5):
    '''
    '''
    rps, Fs = np.ascontiguousarray(rps), np.ascontiguousarray(Fs)
    assert rps.size == Fs.size

    # isolate different regimes
    Fs = Fs*1367*1e7*1e-4   # erg/s/cm^2
    rocky = rps < 1.5
    small = (1.5 <= rps) & (rps < 4)
    neptune = (4 <= rps) & (rps < 13.668)
    giant = rps >= 13.668

    # compute mean mass in each regime
    mps = np.zeros(rps.size)
    mps[rocky]   = .44*rps[rocky]**3 + .614*rps[rocky]**4
    mps[small]   = 2.69 * rps[small]**(.93)
    mps[neptune] = (rps[neptune]*Fs[neptune]**.03 / 1.78)**(1/.53)
    mps[giant]   = np.random.uniform(150,2e3)

    return mps


def _get_sigK(Kdetsig, P, Ms, mp):
    '''
    Compute the desired semi-ampliutde detection measurement uncertainty.
    '''
    K = rvs.RV_K(P, Ms, mp)
    return K / float(Kdetsig)


def _write_results2file(output_fname, magiclistofstuff2write):
    '''
    Write the resulting parameters to a .dat file.
    '''
    # create header with indices
    maglabels = ''
    for i in range(magiclistofstuff2write[3].size):
        maglabels += '%s magnitude\n'%magiclistofstuff2write[9][i]
    hdr = 'Orbital period (days)\nPlanetary radius (Earth radii)\nPlanetary mass(Earth masses)\nRV semi-amplitude(m/s)\n%sStellar mass (Solar masses)\nStellar Radius (Solar radii)\nEffective temperature (K)\n[Fe/H] (Solar units)\nProjected rotation velocity (km/s)\nRotation period (days)\nSpectral resolution\nAperture (meters)\nThroughput\nRV noise floor (m/s)\nReference wavelength (microns)\nTarget SNR\nMaximum fractional telluric absorption\nMinimum exposure time (minutes)\nMaximum exposure time (minutes)\nOverhead (minutes)\nPhoton noise limited RV (m/s)\nRV activity rms (m/s)\nAdditional planet RV rms (m/s)\nEffective RV rms (m/s)\nTarget K measurement uncertainty (m/s)\nExposure time (minutes)\nNumber of RV measurements\nNumber of RV measurements with GP\nTotal observing time (hours)\nTotal observing time with GP (hours)'%maglabels
    hdr, hdrv2 = hdr.split('\n'), ''
    for i in range(len(hdr)):
        hdrv2 += '# %i: %s\n'%(i, hdr[i])

    g = hdrv2
    for i in range(len(magiclistofstuff2write)):
        if i == 3:
            for j in range(magiclistofstuff2write[i].size):
                g += '%.4e\n'%magiclistofstuff2write[i][j]
        elif i == 9:
            pass
        else:
            g += '%.4e\n'%magiclistofstuff2write[i]

    # write .dat file
    try:
        os.mkdir('Results')
    except OSError:
        pass
    f = open('Results/%s.dat'%output_fname, 'w')
    f.write(g)
    f.close()
