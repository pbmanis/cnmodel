from cnmodel.util import expfitting
import numpy as np


def test_fit1():
    fit = expfitting.ExpFitting(nexp=1)
    x = np.linspace(0., 50, 500)
    p = [-50., 4., 5.]
    y = p[0] + p[1]*np.exp(-x/p[2])
    res = fit.fit(x, y, fit.fitpars)
    pr = [float(res[k].value) for k in res.keys()]
    print( '\noriginal: ', p)
    print ('fit res:  ', pr)
    for i, v in enumerate(p):
        assert np.allclose(v, pr[i])

def test_fit2():
    fit = expfitting.ExpFitting(nexp=2)
    x = np.linspace(0., 50, 500)
    p = [-50., 4., 5., 1., 4.5]   # last term is ratio of the two time constants (t2 = delta*t1)
    y = p[0] + p[1]*np.exp(-x/p[2]) + p[3]*np.exp(-x/(p[2]*p[4]))
    res = fit.fit(x, y, fit.fitpars)
    pr = [float(res[k].value) for k in res.keys()]
    print ('\noriginal: ', p)
    print ('fit res:  ', pr)
    # we can only do this approximately for 2 exp fits
    for i, v in enumerate(p): # test each one individually
        assert np.allclose(pr[i]/v, 1.0, atol=1e-4, rtol=1e-2)

if __name__ == '__main__':
    test_fit1()
    test_fit2()