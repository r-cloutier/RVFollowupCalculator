from imports import *
from PyAstronomy.pyasl import foldAt
from compute_sigK_QPGP import compute_sigmaK_GP


def compute_nRV_GP(GPtheta, keptheta, sigRV_phot, sigK_target, duration=100):
    assert len(GPtheta) == 5
    assert len(keptheta) == 2
    a, l, G, Pgp, s = GPtheta
    P, K = keptheta

    # search for target sigma by iterating over Nrv coarsely
    Nrvs = np.arange(10, 1001, 90)
    sigKs = np.zeros(Nrvs.size)
    for i in range(Nrvs.size):
        # get rv activity model
        gp = george.GP(a*(george.kernels.ExpSquaredKernel(l) + \
                          george.kernels.ExpSine2Kernel(G,Pgp)))
        t = _uniform_window_function(duration, Nrvs[i])
        erv = np.repeat(sigRV_phot, t.size)
        gp.compute(t, np.sqrt(erv**2 + s**2))
        rv_act = gp.sample(t)
        ##rv_act *= a / rv_act.std()

        # get planet model
        rv_kep = -K*np.sin(2*np.pi*foldAt(t, P))

        # get total rv signal with noise
        rv = rv_act + rv_kep + np.random.randn(t.size) * sigRV_phot

        # compute sigK
        theta = P, 0, K, a, l, G, Pgp, s
        sigKs[i] = compute_sigmaK_GP(theta, t, rv, erv)

    # fit powerlaw in log space to estimate Nrv for the target sigK
    g = np.isfinite(sigKs)
    p = np.poly1d(np.polyfit(np.log(sigKs[g]), np.log(Nrvs[g]), 1))
    Nrv = np.exp(p(np.log(sigK_target)))
    return Nrv


def _uniform_window_function(duration, Nrv):
    return np.linspace(0, duration, Nrv)


def powerlawfunc(xarr, A, alpha, c):
    xarr = np.ascontiguousarray(xarr)
    return A * xarr.astype('float')**alpha + c
