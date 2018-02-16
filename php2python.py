import os, sys
from RVFollowupCalculator import *


def update_input_files(bands, R, aperture, throughput, floor, wlcen, targetsnr,
                       maxtelluric, mintexp, maxtexp, overhead, sigRVphot,
                       sigRVact, sigRVplanets, sigRVeff, P, rp, mp, mags, Ms,
                       Rs, Teff, Z, vsini, Prot):
    # update spectrograph
    f = open('InputFiles/user_spectrograph_template.in', 'r')
    g = f.read()
    f.close()
    g = g.replace('<<bands>>', ''.join(bands))
    g = g.replace('<<R>>', '%i'%R)
    g = g.replace('<<aperture>>', '%.2f'%aperture)
    g = g.replace('<<throughput>>', '%.2f'%throughput)
    g = g.replace('<<floor>>', '%.2f'%floor)
    g = g.replace('<<wlcen>>', '%.2f'%wlcen)
    g = g.replace('<<targetsnr>>', '%.2f'%targetsnr)
    g = g.replace('<<maxtelluric>>', '%.2f'%maxtelluric)
    g = g.replace('<<mintexp>>', '%.2f'%mintexp)
    g = g.replace('<<maxtexp>>', '%.2f'%maxtexp)
    g = g.replace('<<overhead>>', '%.2f'%overhead)
    h = open('InputFiles/user_spectrograph.in', 'w')
    h.write(g)
    h.close()
    
    # update sigRV
    f = open('user_sigRV_template.in', 'r')
    g = f.read()
    f.close()
    g = g.replace('<<sigRVphot>>', '%.4f'%sigRVphot)
    g = g.replace('<<sigRVact>>', '%.4f'%sigRVact)
    g = g.replace('<<sigRVplanets>>', '%.4f'%sigRVplanets)
    g = g.replace('<<sigRVeff>>', '%.4f'%sigRVeff)
    h = open('InputFiles/user_sigRV.in', 'w')
    h.write(g)
    h.close()

    # update planet
    f = open('user_planet_template.in', 'r')
    g = f.read()
    f.close()
    g = g.replace('<<P>>', '%.4f'%P)
    g = g.replace('<<rp>>', '%.3f'%rp)
    g = g.replace('<<mp>>', '%.3f'%mp)
    h = open('InputFiles/user_planet.in', 'w')
    h.write(g)
    h.close()

    # update star
    f = open('user_star_template.in', 'r')
    g = f.read()
    f.close()
    g = g.replace('<<mags>>', ','.join(['%.2f'%i for i in mags]))
    g = g.replace('<<Ms>>', '%.3f'%Ms)
    g = g.replace('<<Rs>>', '%.3f'%Rs)
    g = g.replace('<<Teff>>', '%i'%Teff)
    g = g.replace('<<Z>>', '%.2f'%Z)
    g = g.replace('<<vsini>>', '%.2f'%vsini)
    g = g.replace('<<Prot>>', '%.2f'%Prot)    
    h = open('InputFiles/user_star.in', 'w')
    h.write(g)
    h.close()


def run_calculator(Kdetsig, Ntrials=1, runGP=True):
    nRV_calculator(Kdetsig,
                   input_planet_fname='user_planet.in',`
                   input_star_fname='user_star.in',
                   input_spectrograph_fname='user_spectrograph.in',
                   input_sigRV_fname='user_sigRV.in',
                   output_fname='RVFCoutput',
                   Ntrials=Ntrials, runGP=runGP)

    
if __name__ == '__main__':
    bands = list(sys.argv[1])
    R = float(sys.argv[2])
    aperture = float(sys.argv[3])
    throughput = float(sys.argv[4])
    floor = float(sys.argv[5])
    wlcen = float(sys.argv[6])
    targetsnr = float(sys.argv[7])
    maxtelluric = float(sys.argv[8])
    mintexp = float(sys.argv[9])
    maxtexp = float(sys.argv[10])
    overhead = float(sys.argv[11])
    sigRVphot = float(sys.argv[12])
    sigRVact = float(sys.argv[13])
    sigRVplanets = float(sys.argv[14])
    sigRVeff = float(sys.argv[15])
    P = float(sys.argv[16])
    rp = float(sys.argv[17])
    mp = float(sys.argv[18])
    mags = list(sys.argv[19])
    Ms = float(sys.argv[20])
    Rs = float(sys.argv[21])
    Teff = float(sys.argv[22])
    Z = float(sys.argv[23])
    vsini = float(sys.argv[24])
    Prot = float(sys.argv[25])
    Kdetsig = float(sys.argv[26])
    Ntrials = float(sys.argv[27])
    runGP = int(sys.argv[28])
    
    update_input_files(bands, R, aperture, throughput, floor, wlcen, targetsnr,
                       maxtelluric, mintexp, maxtexp, overhead, sigRVphot,
                       sigRVact, sigRVplanets, sigRVeff, P, rp, mp, mags, Ms,
                       Rs, Teff, Z, vsini, Prot):

    run_calulator(Kdetsig, Ntrials, runGP)
