from imports import *


def create_pdf(output_fname, magiclistofstuff2write):

    # make sure the output names are not common with any scripts
    assert output_fname not in [s[:-3] for s in glob.glob('*py')]

    # add parameters to tex file
    g = write_pdf_str(output_fname, magiclistofstuff2write)

    # create plots if necessary
    sigRV_acts, sigRV_planets, sigRV_effs,_,_,nRVs, nRVGPs, tobss, tobsGPs = \
                                                    magiclistofstuff2write[-9:]
    if sigRV_acts.size > 1:
        plot_hist(np.array([sigRV_acts, sigRV_planets, sigRV_effs]).T,
                  ['RV activity', 'RV planets', 'Effective RV'],
                  '%s_sigRVs.png'%output_fname)
        plot_hist(np.array([nRVs, nRVGPs]).T, ['nRV (white)', 'nRV (w/ GP)'],
                  '%s_nRVs.png'%output_fname)
        plot_hist(np.array([tobss, tobsGPs]).T,
                  ['Observing time (white)','Observing time (w/ GP)'],
                  '%s_tobs.png'%output_fname)
        figures = True        
        
    # write latex file
    tex_fname = '%s.tex'%output_fname
    h = open(tex_fname, 'w')
    h.write(g)
    h.close()

    # create the pdf and clean up
    os.system('pdflatex %s'%tex_fname)
    os.system('mv %s.pdf Results/'%output_fname)
    os.system('rm %s*'%output_fname)


    
def write_pdf_str(output_fname, magiclistofstuff2write):
    P, rp, mp, K, mags, Ms, Rs, Teff, Z, vsini, Prot, band_strs, R, aperture, throughput, RVnoisefloor, centralwl_microns, SNRtarget, transmission_threshold, texpmin, texpmax, toverhead, sigRV_phot, sigRV_acts, sigRV_planets, sigRV_effs, sigK_target, texp, nRVs, nRVGPs, tobss, tobsGPs = magiclistofstuff2write
    Kdetsig = K / sigK_target
    
    # open/modify latex template
    f = open('outputpdf_template.tex', 'r')
    g = f.read()
    f.close()

    # add parameters
    g = g.replace('<<title>>', output_fname)
    magstr = ', '.join(['%s = %.2f'%(band_strs[i], mags[i])
                        for i in range(len(mags))])
    g = g.replace('<<mags>>', magstr)
    fastlist = ['P', 'rp', 'mp', 'K', 'Ms', 'Rs', 'Teff', 'Z', 'vsini', 'Prot',
                'R', 'aperture', 'throughput', 'RVnoisefloor',
                'centralwl_microns', 'SNRtarget', 'transmission_threshold',
                'texpmin', 'texpmax', 'toverhead', 'sigRV_phot', 'Kdetsig',
                'sigK_target']
    for i in range(len(fastlist)):
        if fastlist[i] in ['Teff','R']:
            fmt = '%i'
        elif fastlist[i] in ['Ms','Rs','sigK_target']:
            fmt = '%.3f'
        elif fastlist[i] in ['Kdetsig']:
            fmt = '%.1f'
        else:
            fmt = '%.2f'
        g = g.replace('<<%s>>'%fastlist[i], fmt%eval(fastlist[i]))

    # add sigRVs depending on size
    if sigRV_acts.size > 1:
        g = g.replace('<<sigRV_act>>',
            'Median RV rms from activity [m/s] & %.2f'%(np.median(sigRV_acts)))
        g = g.replace('<<sigRV_planets>>',
'Median RV rms from additional planets [m/s] & %.2f'%(np.median(sigRV_planets)))
        g = g.replace('<<sigRV_eff>>',
        'Median effective RV precision [m/s] & %.2f'%(np.median(sigRV_effs)))
    else:
        g = g.replace('<<sigRV_act>>',
                      'RV rms from activity [m/s] & %.2f'%float(sigRV_acts))
        g = g.replace('<<sigRV_planets>>',
            'RV rms from additional planets [m/s] & %.2f'%float(sigRV_planets))
        g = g.replace('<<sigRV_eff>>',
        'Effective RV precision [m/s] & %.2f'%float(sigRV_effs))        
        
    return g


def plot_hist(arr, labels, fname):
    Ntrials, Narr = arr.shape
    assert len(labels) == Narr
    
    fig = plt.figure(figsize=(6,4))
    ax = fig.add_subplot(111)
    cols = ['b', 'g', 'r', 'k', 'c']
    for i in range(Narr):
        ax.hist(arr[:,i], bins=int(Ntrials/20), histtype='step',
                color=cols[i], lw=3, label=labels[i])

    ax.set_xlabel(xlabel), ax.set_ylabel('Number of trials')
    ax.legend()
    plt.savefig(fname)
    plt.close('all')
