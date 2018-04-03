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
    

def X2Teff(X, a0, a1, a2, a3):
    '''general relation between color X and Teff'''
    return a0 + a1*X + a2*X**2 + a3*X**3




# TEFF-COLOUR RELATIONS  
def BV2Teff(B_V):
    params = 9084, −7736, 4781, −1342.9
    return X2Teff(B_V, params*)   # Boyajian+2013

def VR2Teff(V_R):
    params = 9335, −9272, 5579, −1302.5
    return X2Teff(V_R, params*)   # Boyajian+2013

def VI2Teff(V_I):
    params = 9354, −7178, 3226, −518.2
    return X2Teff(V_I, params*)   # Boyajian+2013

def VJ2Teff(V_J):
    params = 9052, -3972, 1039, −101.0
    return X2Teff(V_J, params*)   # Boyajian+2013

def VH2Teff(V_H):
    params = 8958, −3023, 632, −52.9
    return X2Teff(V_H, params*)   # Boyajian+2013

def VK2Teff(V_K):
    params = 8984, −2914, 588, −47.4
    return X2Teff(V_K, params*)   # Boyajian+2013

def RJ2Teff(R_J):
    params = 9019, −5767, 2209, −310.3
    return X2Teff(R_J, params)   # Boyajian+2013

def _RH2Teff(R_H):
    params = 9695, −4791, 1432, −175.0
    return X2Teff(R_H, params*)

def _RK2Teff(R_K):
    params = 9683, −4479, 1268, −147.8
    return X2Teff(R_K, params*)





# CONVERSION FUNCTIONS
def V2B(Vmag, Teff):
    B_Vs = np.linspace(-.02, 1.73, 100)
    fint = interp1d(BV2Teff(B_Vs), B_Vs)
    B_V = fint(Teff)
    Bmag = B_V + Vmag
    return Bmag

def V2R(Vmag, Teff):
    V_Rs = np.linspace(0., 1.69, 100)
    fint = interp1d(VR2Teff(V_Rs), V_Rs)
    V_R = fint(Teff)
    Rmag = Vmag - V_R
    return Rmag

def V2I(Vmag, Teff):
    V_Is = np.linspace(-.02., 2.77, 100)
    fint = interp1d(VI2Teff(V_Is), V_Is)
    V_I = fint(Teff)
    Imag = Vmag - V_I
    return Imag

def V2J(Vmag, Teff):
    V_Js = np.linspace(-.12, 4.24, 100)
    fint = interp1d(VJ2Teff(V_Js), V_Js)
    V_J = fint(Teff)
    Jmag = Vmag - V_J
    return Jmag

def V2H(Vmag, Teff):
    V_Hs = np.linspace(-.13, 4.77, 100)
    fint = interp1d(VJ2Teff(V_Hs), V_Hs)
    V_H = fint(Teff)
    Hmag = Vmag - V_H
    return Hmag

def V2K(Vmag, Teff):
    V_Ks = np.linspace(-.15, 5.04, 100)
    fint = interp1d(VJ2Teff(V_Ks), V_Ks)
    V_K = fint(Teff)
    Kmag = Vmag - V_K
    return Kmag

def J2H(Jmag, Teff):
    R_Js = np.linspace(.09, 2.58, 100)
    fint = interp1d(RJ2Teff(R_Js), R_Js)
    R_J = fint(Teff)
    Rmag = R_J + Jmag
    R_Hs = np.linspace(.07, 3.17, 100)
    fint = interp1d(RH2Teff(R_Hs), R_Hs)
    R_H = fint(Teff)
    Hmag = Rmag - R_H
    return Hmag

def J2K(Jmag, Teff):
    R_Js = np.linspace(.09, 2.58, 100)
    fint = interp1d(RJ2Teff(R_Js), R_Js)
    R_J = fint(Teff)
    Rmag = R_J + Jmag
    R_Ks = np.linspace(.06, 3.43, 100)
    fint = interp1d(RK2Teff(R_Ks), R_Ks)
    R_K = fint(Teff)
    Kmag = Rmag - R_K
    return Kmag
