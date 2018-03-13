import numpy as np
import spiceypy as spice
import pandas as pd
import scipy as sp
import math as mt
from scipy.optimize.nonlin import NoConvergence 
from scipy.interpolate import interp1d 
import matplotlib.pyplot as plt
from scipy import signal as sg

print('Start of Program.\n')

def nonrelbroyfunc(nonrelguessx): 
    global passdf,passc,passf,passbetaa,passvra,passvza,passdeltab,passvrb,passvzb,passzb,passra,passza,passgamma
    return [passdf * passc / passf - np.cos(passbetaa) * passvra - np.sin(passbetaa) * passvza + np.sin(passdeltab) * passvrb + np.cos(passdeltab) * passvzb + np.cos(passbetaa - nonrelguessx[1]) * passvra + np.sin(passbetaa - nonrelguessx[1]) * passvza - np.sin(passdeltab - nonrelguessx[0]) * passvrb - np.cos(passdeltab - nonrelguessx[0]) * passvzb, (-1) * passzb * np.sin(passdeltab - nonrelguessx[0]) - np.sqrt(passra**2 + passza**2) * np.sin(passbetaa - nonrelguessx[1] - passgamma)]

def relbroyfunc(relguessx): 
    global passdf,passc,passf,passbetaa,passvra,passvza,passdeltab,passvrb,passvzb,passzb,passra,passza,passgamma
    return [passdf / passf - (1 + passvrb * np.sin(passdeltab-relguessx[0]) / passc + passvzb * np.cos(passdeltab-relguessx[0]) / passc + 0.5 * passvbsqd / passc**2 - passub/passc**2) / (1 + passvra * np.cos(passbetaa-relguessx[1]) / passc + passvza * np.sin(passbetaa-relguessx[1]) / passc + 0.5 * passvasqd / passc**2 - passua/passc**2) + (1 + passvrb * np.sin(passdeltab) / passc + passvzb * np.cos(passdeltab) / passc + 0.5 * passvbsqd / passc**2 - passub/passc**2) / (1 + passvra * np.cos(passbetaa) / passc + passvza * np.sin(passbetaa) / passc + 0.5 * passvasqd / passc**2 - passua/passc**2),(-1) * passzb * np.sin(passdeltab - relguessx[0]) - np.sqrt(passra**2 + passza**2) * np.sin(passbetaa - relguessx[1] - passgamma)]

nconvergence =  1000
freq         =  8423*pow(10,6)

for iobs in range(0,nobs):
    print(iobs)
    passdf = deltaf[iobs]
    passf = freq
    passc = c_mks
    passbetaa = betaa
    passdeltab = deltab
    passvra = vra
    passvza = vza
    passvrb = vrb
    passvzb = vzb
    passzb = zb
    passra = ra
    passza = za
    passgamma = gamma

    nonrelguessx = [0,0] 
    try :
        nonrelresult = sp.optimize.broyden1( nonrelbroyfunc, nonrelguessx, maxiter = 5000, x_tol = 10**-20)
        converged = True
    except NoConvergence as e:
        nonrelresult = e.args[0]
        converged = False

    nonreldeltararr[iobs] = nonrelresult[0]
    nonrelbetararr[iobs] = nonrelresult[1]

    nonrelguessx = [0,0] 
    try :
        relresult = sp.optimize.broyden1( relbroyfunc, relguessx, maxiter = 5000, x_tol = 10**-20)
        converged = True
    except NoConvergence as e:
        relresult = e.args[0]
        converged = False

    reldeltararr[iobs] = relresult[0]
    relbetararr[iobs] = relresult[1]
