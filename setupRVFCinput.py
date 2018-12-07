from imports import *

global basedir, fname_star, fname_planet, fname_spec, fname_sigRV
basedir = 'InputFiles'
fname_star = '%s/user_star_template.in'%basedir
fname_planet = '%s/user_planet_template.in'%basedir
fname_spec = '%s/user_spectrograph_template.in'%basedir
fname_sigRV = '%s/user_sigRV_template.in'%basedir


def setup_star(prefix, mag, Ms, Rs, Teff, Z, vsini, Prot):
    '''Setup the RVFC stellar input file given input stellar parameters.'''
    # get template file
    f = open(fname_star, 'r')
    g = f.read()
    f.close()

    # update input stellar parameters
    g = g.replace('<<mag>>', '%.3f'%mag)
    g = g.replace('<<Ms>>', '%.3f'%Ms)
    g = g.replace('<<Rs>>', '%.3f'%Rs)
    g = g.replace('<<Teff>>', '%.3f'%Teff)
    g = g.replace('<<Z>>', '%.3f'%Z)
    g = g.replace('<<vsini>>', '%.3f'%vsini)
    g = g.replace('<<Prot>>', '%.3f'%Prot)

    # write new file
    fname_out = '%s/%s_star.in'%(basedir, prefix)
    h = open(fname_out, 'w')
    h.write(g)
    h.close()

    return fname_out



def setup_planet(prefix, P, rp, mp=None):
    '''Setup the RVFC planet input file given input planet parameters.'''
    # get template file
    f = open(fname_planet, 'r')
    g = f.read()
    f.close()

    # update input planet parameters
    g = g.replace('<<P>>', '%.8f'%P)
    g = g.replace('<<rp>>', '%.3f'%rp)
    mp = float(mp) if mp != None else 0.
    g = g.replace('<<mp>>', '%.3f'%mp)

    # write new file
    fname_out = '%s/%s_planet.in'%(basedir, prefix)
    h = open(fname_out, 'w')
    h.write(g)
    h.close()

    return fname_out



def setup_spectrograph(prefix, wlmin, wlmax, R, aperture, throughput, RVfloor,
                       maxtelluric, toverhead):
    '''Setup the RVFC spectrograph input file given input parameters.'''
    # get template file
    f = open(fname_spec, 'r')
    g = f.read()
    f.close()

    # update input spectrograph parameters
    g = g.replace('<<wlmin>>', '%.1f'%wlmin)
    g = g.replace('<<wlmax>>', '%.1f'%wlmax)
    g = g.replace('<<R>>', '%i'%R)
    g = g.replace('<<aperture>>', '%.3f'%aperture)
    g = g.replace('<<throughput>>', '%.4f'%throughput)
    g = g.replace('<<RVnoisefloor>>', '%.3f'%RVfloor)
    g = g.replace('<<maxtelluric>>', '%.4f'%maxtelluric)
    g = g.replace('<<toverhead>>', '%.2f'%toverhead)

    # write new file
    fname_out = '%s/%s_spectrograph.in'%(basedir, prefix)

    h = open(fname_out, 'w')
    h.write(g)
    h.close()

    return fname_out



def setup_sigRV(prefix, texp_sec, sigphot=None, sigact=None, sigplanets=None,
                sigeff=None):
    '''Setup the RVFC RV noise input file given input parameters.'''
    # get template file
    f = open(fname_sigRV, 'r')
    g = f.read()
    f.close()

    # update input RV noise parameters
    g = g.replace('<<texp>>', '%.3f'%(texp_sec/60.))
    sigphot = float(sigphot) if sigphot != None else 0.
    g = g.replace('<<sigphot>>', '%.3f'%sigphot)
    sigact = float(sigact) if sigact != None else -1.
    g = g.replace('<<sigact>>', '%.3f'%sigact)
    sigplanets = float(sigplanets) if sigplanets != None else -1.
    g = g.replace('<<sigplanets>>', '%.3f'%sigplanets)
    sigeff = float(sigeff) if sigeff != None else 0.
    g = g.replace('<<sigeff>>', '%.3f'%sigeff)

    # write new file
    fname_out = '%s/%s_sigRV.in'%(basedir, prefix)
    h = open(fname_out, 'w')
    h.write(g)
    h.close()

    return fname_out


def clean_files(prefix):
    '''Remove the input files created by any of the above functions.'''
    fs = glob.glob('%s/%s_*in'%(basedir, prefix))
    for f in fs:
	os.system('rm %s'%f)
