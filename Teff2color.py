from imports import *


# MISC FUNCTIONS
def define_band_dictionaries():
    # ergs/s/cm2/micron
    zeropoints = {'U': 3.678e-5, 'B': 6.293e-5,
                  'V': 3.575e-5, 'R': 1.882e-5,
                  'I': 9.329e-6, 'Y': 5.949e-6,
                  'J': 2.985e-6, 'H': 1.199e-6,
                  'K': 4.442e-7}
    # nm
    wlcens = {'U': 353.1, 'B': 443.0,
              'V': 553.7, 'R': 694.0,
              'I': 878.1, 'Y': 1025.9,
              'J': 1254.5, 'H': 1631.0,
              'K': 2149.8}
    return wlcens, zeropoints
    

def X2Teff_AFGK(X, *args):
    '''general relation between color X and Teff for AFGK stars from
    2013ApJ...771...40B'''
    a0, a1, a2, a3 = args 
    return a0 + a1*X + a2*X**2 + a3*X**3


def X2Teff_M(X, FeH, *args):
    '''general relation between color X and Teff for M dwarfs from 
    2012ApJ...757..112B'''
    a0, a1, a2, a3, a4, a5 = args
    return a0 + a1*X + a2*X**2 + a3*X*FeH + a4*FeH + a5*FeH**2



# TEFF-COLOUR RELATIONS  
def BV2Teff(B_V, FeH=None):
    params_AFGK = 9084, -7736, 4781, -1342.9
    params_M = 8010, -4095, 819, 133, 39, -362
    return X2Teff_AFGK(B_V, *params_AFGK) if FeH is None \
        else X2Teff_M(B_V, FeH, *params_M)

def VR2Teff(V_R, FeH=None):
    params_AFGK = 9335, -9272, 5579, -1302.5
    params_M = 7646, -4295, 1058, 304, -77, 102
    return X2Teff_AFGK(V_R, *params_AFGK) if FeH is None \
        else X2Teff_M(V_R, FeH, *params_M)
    
def VI2Teff(V_I, FeH=None):
    params_AFGK = 9354, -7178, 3226, -518.2
    params_M = 7325, -2262, 313, -15, 393, 733
    return X2Teff_AFGK(V_I, *params_AFGK) if FeH is None \
        else X2Teff_M(V_I, FeH, *params_M)

def VJ2Teff(V_J, FeH=None):
    params_AFGK = 9052, -3972, 1039, -101.0
    params_M = 7308, -1775, 198, 71, 100, 317
    return X2Teff_AFGK(V_J, *params_AFGK) if FeH is None \
        else X2Teff_M(V_J, FeH, *params_M)
    
def VH2Teff(V_H, FeH=None):
    params_AFGK = 8958, -3023, 632, -52.9
    params_M = 7641, -1611, 151, 177, -319, 185
    return X2Teff_AFGK(V_H, *params_AFGK) if FeH is None \
        else X2Teff_M(V_H, FeH, *params_M)

def VK2Teff(V_K, FeH=None):
    params_AFGK = 8984, -2914, 588, -47.4
    params_M = 7643, -1523, 134, 137, -202, 157
    return X2Teff_AFGK(V_K, *params_AFGK) if FeH is None \
        else X2Teff_M(V_K, FeH, *params_M)

#def RJ2Teff(R_J, FeH=None):
#    params = 9019, -5767, 2209, -310.3
#    return X2Teff(R_J, params)   # Boyajian+2013
#def RH2Teff(R_H, FeH=None):
#    params = 9695, -4791, 1432, -175.0
#    return X2Teff(R_H, params*)
#def RK2Teff(R_K, FeH=None):
#    params = 9683, -4479, 1268, -147.8
#    return X2Teff(R_K, params*)




# CONVERSION FUNCTIONS
def V2B(Vmag, Teff, FeH):
    B_Vs = np.linspace(-.02, 1.73, 100)
    FeH = float(FeH) if Teff < 4e3 else None
    fint = interp1d(BV2Teff(B_Vs, FeH), B_Vs)
    B_V = fint(Teff)
    Bmag = B_V + Vmag
    return Bmag

def V2R(Vmag, Teff, FeH):
    V_Rs = np.linspace(0., 1.69, 100)
    FeH = float(FeH) if Teff < 4e3 else None
    fint = interp1d(VR2Teff(V_Rs, FeH), V_Rs)
    V_R = fint(Teff)
    Rmag = Vmag - V_R
    return Rmag

def V2I(Vmag, Teff, FeH):
    V_Is = np.linspace(-.02, 2.77, 100)
    FeH = float(FeH) if Teff < 4e3 else None
    fint = interp1d(VI2Teff(V_Is, FeH), V_Is)
    V_I = fint(Teff)
    Imag = Vmag - V_I
    return Imag

def V2J(Vmag, Teff, FeH):
    V_Js = np.linspace(-.12, 4.24, 100)
    FeH = float(FeH) if Teff < 4e3 else None
    fint = interp1d(VJ2Teff(V_Js, FeH), V_Js)
    V_J = fint(Teff)
    Jmag = Vmag - V_J
    return Jmag

def V2H(Vmag, Teff, FeH):
    V_Hs = np.linspace(-.13, 4.77, 100)
    FeH = float(FeH) if Teff < 4e3 else None
    fint = interp1d(VJ2Teff(V_Hs, FeH), V_Hs)
    V_H = fint(Teff)
    Hmag = Vmag - V_H
    return Hmag

def V2K(Vmag, Teff, FeH):
    V_Ks = np.linspace(-.15, 5.04, 100)
    FeH = float(FeH) if Teff < 4e3 else None
    fint = interp1d(VJ2Teff(V_Ks, FeH), V_Ks)
    V_K = fint(Teff)
    Kmag = Vmag - V_K
    return Kmag

def J2V(Jmag, Teff, FeH):
    V_Js = np.linspace(-.12, 4.24, 100)
    FeH = float(FeH) if Teff < 4e3 else None
    fint = interp1d(VJ2Teff(V_Js, FeH), V_Js)
    V_J = fint(Teff)
    Vmag = V_J + Jmag
    return Vmag    

def J2H(Jmag, Teff, FeH):
    V_Js = np.linspace(-.12, 5, 100)
    FeH = float(FeH) if Teff < 4e3 else None
    assert Teff >= VJ2Teff(V_Js, FeH).min()
    fint = interp1d(VJ2Teff(V_Js, FeH), V_Js)
    assert
    V_J = fint(Teff)
    Vmag = V_J + Jmag
    Hmag =  V2H(Vmag, Teff, FeH)
    return Hmag

def J2K(Jmag, Teff, FeH):
    V_Js = np.linspace(-.12, 4.24, 100)
    FeH = float(FeH) if Teff < 4e3 else None
    fint = interp1d(VJ2Teff(V_Js, FeH), V_Js)
    V_J = fint(Teff)
    Vmag = V_J + Jmag
    Kmag =  V2K(Vmag, Teff, FeH)
    return Kmag


def V2all(Vmag, Teff, FeH):
    # return all magnitudes
    Umag = np.nan
    Bmag = V2B(Vmag, Teff, FeH)
    Rmag = V2R(Vmag, Teff, FeH)
    Imag = V2I(Vmag, Teff, FeH)
    Jmag = V2J(Vmag, Teff, FeH)
    Hmag = V2H(Vmag, Teff, FeH)
    Kmag = V2K(Vmag, Teff, FeH)
    Ymag = np.nan #np.mean([Imag,Jmag])
    return Umag, Bmag, Vmag, Rmag, Imag, Ymag, Jmag, Hmag, Kmag

def J2all(Jmag, Teff, FeH):
    # return all magnitudes
    Umag = np.nan
    Bmag = V2B(J2V(Jmag, Teff, FeH), Teff, FeH)
    Vmag = J2V(Jmag, Teff, FeH)
    Rmag = V2R(Vmag, Teff, FeH)
    Imag = V2I(Vmag, Teff, FeH)
    Ymag = np.nan
    Hmag = J2H(Jmag, Teff, FeH)
    Kmag = J2K(Jmag, Teff, FeH)
    return Umag, Bmag, Vmag, Rmag, Imag, Ymag, Jmag, Hmag, Kmag
