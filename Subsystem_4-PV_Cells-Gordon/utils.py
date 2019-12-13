############################
# Created on Wed Dec  4 15:51:06 2019
# author: Gordon Cheung Yat Hei
# CID: 01083012
# Project: Optimisation
############################

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize
import sympy as sym


## All the fitting functions defined to make easy usage

def fit_sin(tt, yy):
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

def fit_log(tt, yy):
    tt = np.array(tt)
    yy = np.array(yy)
    def log_func(x, a, b, c):
        return a * np.log(b * x) + c

    popt, pcov = optimize.curve_fit(log_func, tt, yy)
    a, b, c = popt
    fitfunc = lambda x: a * np.log(b * x) + c
    return {"a": a, "b": b, "c": c,
            "fitfunc": fitfunc, "maxcov": np.max(pcov), "rawres": (popt,pcov)}

def fit_exp(tt, yy):
    tt = np.array(tt)
    yy = np.array(yy)
    def exp_func(x, a, b, c):
        return a * np.exp(b * x) + c
    popt, pcov = optimize.curve_fit(exp_func, tt, yy)
    a, b, c = popt
    fitfunc = lambda x: a * np.exp(b * x) + c
    return {"a": a, "b": b, "c": c,
            "fitfunc": fitfunc, "maxcov": np.max(pcov), "rawres": (popt,pcov)}
    
def fit_lin(tt, yy):
    tt = np.array(tt)
    yy = np.array(yy)
    def exp_func(x, a, b):
        return a * x + b
    popt, pcov = optimize.curve_fit(exp_func, tt, yy)
    a, b = popt
    fitfunc = lambda x: a * x + b
    return {"a": a, "b": b,
            "fitfunc": fitfunc, "maxcov": np.max(pcov), "rawres": (popt,pcov)}




def roof_width(width, beta):
    # compute the width of the roof using the length and the angle
    return (width*np.sin(beta))/np.sin(np.pi-2*beta)

def max_tiling(p_length, p_width, r_length, r_width):
    # function to compute that maximum amount of panels that can be put on the roof
    # assuming all panels to be placed need to be the same model
    if r_width >= p_length:
        # if the panels can be fitted perpendicularly
        number = np.floor(r_width/p_length)*np.floor(r_length/p_width)
        if r_width-p_length >= p_width:
            # if more panels can be fitted parallely
            number = number + np.floor(r_width-p_length/p_width)*np.floor(r_length/p_length)
    elif r_width >= p_width:
        # otherwise if we can fit the panels parallely
        number = np.floor(r_width/p_width)*np.floor(r_length/p_length)
    else: number = 0
    return number

def bills(energy, consumption, domain):
    # function to compute the     
    buy_price = 12.5/100 #grid power selling price/kWh
#    sell_price = 5.5/100 #grid power buying price/kWh
    sell_price = 12.5/100 #grid power buying price/kWh
    cost = 0 # initialising costs
    sell = 0
    for i in np.arange(domain[0],domain[1]):
        difference = consumption(i)-energy(i)
        if difference > 0:
            # if there is more consumption than energy generated
            cost = cost + difference*buy_price
        else:
            # if the excess energy can be sold back to the grid
            sell = sell + difference*sell_price
    return cost-sell   

def pareto_set(name, data_1, data_2, min1=True, min2=True):
    # function to find the pareto set in the set of data points
    # limited to 2 dimensions, can be further expanded
    if min1 == True:
        index = data_1.argsort()[::1]
        data_1 = data_1[index] #sort from largest to smallest
        data_2 = data_2[index]
        name = name[index]
    else: 
        index = data_1.argsort()[::-1]
        data_1 = data_1[index] #sort from smallest to largest
        data_2 = data_2[index]
        name = name[index[::1]]
    pareto = [index[len(index)-1]]
    name_set = [name[index[0]]]
    ref = data_2[index[0]]
    for i in np.arange(1,len(data_1)):
        if min2==True:
            if data_2[index[i]]<ref:
                pareto.append(index[len(index)-1-i])
                ref = data_2[index[i]]
                name_set.append(name[index[i]])
        else:
            if data_2[index[i]]>ref:
                pareto.append(index[len(index)-1-i])
                ref = data_2[index[i]]
                name_set.append(name[index[i]])
    return pareto, name_set
 