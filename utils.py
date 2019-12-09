# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 15:51:06 2019

@author: Gordon Cheung
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize
import sympy as sym

def fit_sin(tt, yy):
    t = sym.Symbol('t')
    # algorithm taken from:
    # https://stackoverflow.com/questions/16716302/how-do-i-fit-a-sine-curve-to-my-data-with-pylab-and-numpy#16716964
    ## Fit sin to the input time sequence, and return fitting parameters "A", "omega", "b", "c" and fitfunc"
    tt = np.array(tt)
    yy = np.array(yy)
    ff = np.fft.fftfreq(len(tt), (tt[1]-tt[0]))   # assume uniform spacing
    Fyy = abs(np.fft.fft(yy))
    guess_freq = abs(ff[np.argmax(Fyy[1:])+1])   # excluding the zero frequency "peak", which is related to offset
    guess_amp = np.std(yy) * 2.**0.5
    guess_offset = np.mean(yy)
    guess = np.array([guess_amp, 2.*np.pi*guess_freq, 0., guess_offset])

    def sinfunc(t, A, w, b, c):  
        return A * np.sin(w*t + b) + c

    popt, pcov = optimize.curve_fit(sinfunc, tt, yy, p0=guess)
    A, w, b, c = popt
    fitfunc = lambda t: A * np.sin(w*t + b) + c
    return {"A": A, "w": w, "b": b, "c": c,
            "fitfunc": fitfunc, "maxcov": np.max(pcov), "rawres": (guess,popt,pcov)}