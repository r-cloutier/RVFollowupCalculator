from imports import *

def get_data():
    '''get theoretical data from Lejeune+1998'''
    d = np.loadtxt('input_data/asu.tsv',delimiter=';', skiprows=51).T
    g = d[0] <= 1e4
    Teff,logg,FeH,U_B,B_V,V_R,V_I,V_K,R_I,J_H,J_K,J_L = d[:,g]
    return Teff, logg, FeH, U_B, B_V, V_R, V_I, V_K, R_I, J_H, J_K, J_L

def isolate_logg_FeH(logg, FeH):
    _,loggs,FeHs,_,_,_,_,_,_,_,_,_ = get_data()
    g = (abs(loggs-logg) == np.min(abs(loggs-logg))) & \
        (abs(FeHs-FeH) == np.min(abs(FeHs-FeH)))
    return np.where(g)[0]


# MAGNITUDE CONVERSIONS
def _V2U(Vmag, Teff, logg, FeH):
    Teffs,_,_,U_Bs,B_Vs,_,_,_,_,_,_,_ = get_data()
    g = isolate_logg_FeH(logg, FeH)
    #plt.scatter(Teffs[g], B_Vs[g]), plt.show()
    Bmag = _V2B(Vmag, Tedd, logg, FeH)
    fint = interp1d(Teffs, U_Bs)
    U_B = fint(Teff)
    Umag = U_B + Bmag
    return Umag

def _V2B(Vmag, Teff, logg, FeH):
    Teffs,_,_,_,B_Vs,_,_,_,_,_,_,_ = get_data()
    g = isolate_logg_FeH(logg, FeH)
    #plt.scatter(Teffs[g], B_Vs[g]), plt.show()
    fint = interp1d(Teffs, B_Vs)
    B_V = fint(Teff)
    Bmag = B_V + Vmag
    return Bmag

def _V2R(Vmag, Teff, logg, FeH):
    Teffs,_,_,_,_,V_Rs,_,_,_,_,_,_ = get_data()
    g = isolate_logg_FeH(logg, FeH)
    #plt.scatter(Teffs[g], V_Rs[g]), plt.show()
    fint = interp1d(Teffs, V_Rs)
    V_R = fint(Teff)
    Rmag = Vmag - V_R
    return Rmag

def _V2I(Vmag, Teff, logg, FeH):
    Teffs,_,_,_,_,_,V_Is,_,_,_,_,_ = get_data()
    g = isolate_logg_FeH(logg, FeH)
    #plt.scatter(Teffs[g], V_Is[g]), plt.show()
    fint = interp1d(Teffs, V_Is)
    V_I = fint(Teff)
    Imag = Vmag - V_I
    return Imag

def _V2J(Vmag, Teff, logg, FeH):
    Kmag = _V2K(Vmag, Teff, logg, FeH)
    Teffs,_,_,_,_,_,_,V_Ks,_,_,J_Ks,_ = get_data()
    g = isolate_logg_FeH(logg, FeH)
    #plt.scatter(Teffs[g], V_Ks[g]), plt.show()
    fint = interp1d(Teffs, J_Ks)
    J_K = fint(Teff)
    Jmag = J_K + Kmag
    return Jmag

def _V2H(Vmag, Teff, logg, FeH):
    Jmag = _V2J(Vmag, Teff, logg, FeH)
    Teffs,_,_,_,_,_,_,_,_,J_Hs,_,_ = get_data()
    g = isolate_logg_FeH(logg, FeH)
    #plt.scatter(Teffs[g], J_Hs[g]), plt.show()
    fint = interp1d(Teffs, J_Hs)
    J_H = fint(Teff)
    Hmag = Jmag - J_H
    return Hmag

def _V2K(Vmag, Teff, logg, FeH):
    Teffs,_,_,_,_,_,_,V_Ks,_,_,_,_ = get_data()
    g = isolate_logg_FeH(logg, FeH)
    #plt.scatter(Teffs[g], V_Ks[g]), plt.show()
    fint = interp1d(Teffs, V_Ks)
    V_K = fint(Teff)
    Kmag = Vmag - V_K
    return Kmag

def _J2U(Jmag, Teff, logg, FeH):
    Teffs,_,_,_,_,_,_,V_Ks,_,_,J_Ks,_ = get_data()
    g = isolate_logg_FeH(logg, FeH)
    #plt.scatter(Teffs[g], V_Ks[g]), plt.show()
    fint = interp1d(Teffs, J_Ks)
    J_K = fint(Teff)
    Kmag = Jmag - J_K
    fint = interp1d(Teffs, V_Ks)
    V_K = fint(Teff)
    Vmag = V_K + Kmag
    Umag = _V2U(Vmag, Teff, logg, FeH)
    return Umag

def _J2B(Jmag, Teff, logg, FeH):
    Teffs,_,_,_,_,_,_,V_Ks,_,_,J_Ks,_ = get_data()
    g = isolate_logg_FeH(logg, FeH)
    #plt.scatter(Teffs[g], V_Ks[g]), plt.show()
    fint = interp1d(Teffs, J_Ks)
    J_K = fint(Teff)
    Kmag = Jmag - J_K
    fint = interp1d(Teffs, V_Ks)
    V_K = fint(Teff)
    Vmag = V_K + Kmag
    Bmag = _V2B(Vmag, Teff, logg, FeH)
    return Bmag

def _J2V(Jmag, Teff, logg, FeH):
    Teffs,_,_,_,_,_,_,V_Ks,_,_,J_Ks,_ = get_data()
    g = isolate_logg_FeH(logg, FeH)
    #plt.scatter(Teffs[g], V_Ks[g]), plt.show()
    fint = interp1d(Teffs, J_Ks)
    J_K = fint(Teff)
    Kmag = Jmag - J_K
    fint = interp1d(Teffs, V_Ks)
    V_K = fint(Teff)
    Vmag = V_K + Kmag
    return Vmag

def _J2R(Jmag, Teff, logg, FeH):
    Teffs,_,_,_,_,_,_,V_Ks,_,_,J_Ks,_ = get_data()
    g = isolate_logg_FeH(logg, FeH)
    #plt.scatter(Teffs[g], V_Ks[g]), plt.show()
    fint = interp1d(Teffs, J_Ks)
    J_K = fint(Teff)
    Kmag = Jmag - J_K
    fint = interp1d(Teffs, V_Ks)
    V_K = fint(Teff)
    Vmag = V_K + Kmag
    Rmag = _V2R(Vmag, Teff, logg, FeH)
    return Rmag

def _J2I(Jmag, Teff, logg, FeH):
    Teffs,_,_,_,_,_,_,V_Ks,_,_,J_Ks,_ = get_data()
    g = isolate_logg_FeH(logg, FeH)
    #plt.scatter(Teffs[g], V_Ks[g]), plt.show()
    fint = interp1d(Teffs, J_Ks)
    J_K = fint(Teff)
    Kmag = Jmag - J_K
    fint = interp1d(Teffs, V_Ks)
    V_K = fint(Teff)
    Vmag = V_K + Kmag
    Imag = _V2I(Vmag, Teff, logg, FeH)
    return Imag

def _J2H(Jmag, Teff, logg, FeH):
    Teffs,_,_,_,_,_,_,_,_,J_Hs,_,_ = get_data()
    g = isolate_logg_FeH(logg, FeH)
    #plt.scatter(Teffs[g], V_Ks[g]), plt.show()
    fint = interp1d(Teffs, J_Hs)
    J_H = fint(Teff)
    Hmag = Jmag - J_H
    return Hmag

def _J2K(Jmag, Teff, logg, FeH):
    Teffs,_,_,_,_,_,_,_,_,_,J_Ks,_ = get_data()
    g = isolate_logg_FeH(logg, FeH)
    #plt.scatter(Teffs[g], V_Ks[g]), plt.show()
    fint = interp1d(Teffs, J_Ks)
    J_K = fint(Teff)
    Kmag = Jmag - J_K
    return Kmag


def V2all(Vmag, Teff, logg, FeH):
    Vmag = float(Vmag)
    Umag = _V2U(Vmag, Teff, logg, FeH)
    Bmag = _V2B(Vmag, Teff, logg, FeH)
    Rmag = _V2R(Vmag, Teff, logg, FeH)
    Imag = _V2I(Vmag, Teff, logg, FeH)
    Ymag = np.nan
    Jmag = _V2J(Vmag, Teff, logg, FeH)
    Hmag = _V2H(Vmag, Teff, logg, FeH)
    Kmag = _V2K(Vmag, Teff, logg, FeH)
    return Umag, Bmag, Vmag, Rmag, Imag, Ymag, Jmag, Hmag, Kmag

def J2all(Jmag, Teff, logg, FeH):
    Jmag = float(Jmag)
    Umag = _J2U(Jmag, Teff, logg, FeH)
    Bmag = _J2B(Jmag, Teff, logg, FeH)
    Vmag = _J2V(Jmag, Teff, logg, FeH)
    Rmag = _J2R(Jmag, Teff, logg, FeH)
    Imag = _J2I(Jmag, Teff, logg, FeH)
    Ymag = np.nan
    Hmag = _J2H(Jmag, Teff, logg, FeH)
    Kmag = _J2K(Jmag, Teff, logg, FeH)
    return Umag, Bmag, Vmag, Rmag, Imag, Ymag, Jmag, Hmag, Kmag
