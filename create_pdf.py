from imports import *

def create_pdf(output_fname, magiclistofstuff2write):
    g = write_pdf_str(output_fname, magiclistofstuff2write)

    # write latex file
    tex_fname = '%s.tex'%output_fname
    h = open(tex_fname, 'w')
    h.write(g)
    h.close()

    # create the pdf
    os.system('pdflatex %s'%tex_fname)
    os

    
def write_pdf_str(output_fname, magiclistofstuff2write):
    P, rp, mp, K, mags, Ms, Rs, Teff, Z, vsini, band_strs, R, aperture, throughput, RVnoisefloor, centralwl_microns, SNRtarget, transmission_threshold, texpmin, texpmax, toverhead, sigRV_phot, sigRV_act, sigRV_planets, sigRV_eff, sigK_target, nRV, nRVGP, texp, tobs, tobsGP = magiclistofstuff2write
    Kdetsig = K / sigK_target
    
    # open/modify latex template
    f = open('outputpdf_template.tex', 'r')
    g = f.read()
    f.close()

    g = g.replace('<<title>>', output_fname)
    magstr = ', '.join(['%s = %.2f'%(band_strs[i], mags[i])
                        for i in range(len(mags))])
    g = g.replace('<<mags>>', magstr)


    return g
