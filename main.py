# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 15:51:29 2019

@author: Gordon Cheung
"""
import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
from solar_irradiation import Solar_irradiation
from utils import *

######## Part 1. Formulation of the problem ########

### Parameters
resolution = 0.1
days = 365 # number of days taken account of in this calculation
latitude = np.deg2rad(50.41)
t = sym.Symbol('t') # time parameter for the optimisation problem

### Part 1.1 Setting up the solar irradiance data
solar = Solar_irradiation(latitude, np.arange(1,days*24, resolution), resolution)

## get the daily averages of solar irradiation
yavg = solar.get_daily_average()

# fit a surrogate model onto the daily average to simplify the algorithm
tt = np.arange(days)
solar_daily = fit_sin(tt, yavg)
#print( "Amplitude=%(amp)s \n Angular freq.=%(omega)s \n phase=%(phase)s, \
#      \n offset=%(offset)s \n Max. Cov.=%(maxcov)s" % solar_daily )
solar_daily_fit = solar_daily["A"]*sym.sin(solar_daily["w"]*t+solar_daily["b"])+solar_daily["c"]
print(sym.diff(solar_daily_fit,t))
plt.figure()
plt.title("Sinusodal curve fitting of the daily solar irradiance")
plt.plot(tt, yavg, "ob", label="daily averages", linewidth=2)
plt.plot(tt, solar_daily["fitfunc"](tt), "y-", label="fit curve", linewidth=2)
plt.legend(loc="best")
plt.show()

## compute the maximum and minimum days for later optimisation
max_day = np.argmax(yavg)
min_day = np.argmin(yavg)
#print(" max day: ", max_day, ", value: ", yavg[max_day])
#print(" min day: ", min_day, ", value: ", yavg[min_day])
yId = solar.get_data()
#plt.figure()
#plt.plot(yId)

# Using the maximum and minimum to slice the days out for later optimisation
max_data = yId[int(max_day*24/resolution):int((max_day*24+23)/resolution)]
min_data = yId[int(min_day*24/resolution):int((min_day*24+23)/resolution)]
#min_data = yId[int(min_day/resolution):int((min_day+23)/resolution)]


## fit a surrogate model onto the daily average to simplify the algorithm
tt = np.arange(0, 23,resolution)
solar_max = fit_sin(tt, max_data)
solar_min = fit_sin(tt, min_data)
#plt.figure()
#plt.title("Sinusodal curve fitting of the day with max and min solar irradiance")
#plt.plot(tt, max_data, "bo", label="max data", linewidth=1)
#plt.plot(tt, solar_max["fitfunc"](tt), "b-", label="max fit", linewidth=2)
#plt.plot(tt, min_data, "ro", label="min data", linewidth=2)
#plt.plot(tt, solar_min["fitfunc"](tt), "r-", label="min fit", linewidth=1)
#plt.legend(loc="best")
#plt.show()


### Part 1.2 Setting up the consumption data 

# Using dataset taken online, with the power consumption adjusted assuming 
# the household in study is average and varies with the national average
# Further adjusted by the region defined using the 2017 consumption data
df = pd.read_csv("C:\\Users\\Gordon Cheung\\OneDrive - Imperial College London\\1.Imperial\\\
                 Year_4\\Optimisation\\Coursework\\Code\\Data\\energy_consumption.csv")
# converting the time scale from days to hours (to be consistent from above)
df['hour']=df['Day']*24
plt.plot(df['hour'],df['Power (adjusted)']) # visualising 

# fit a surrogate model onto the daily average to simplify the algorithm
tt = df['hour']
consumption = fit_sin(df['hour'], df['Power (adjusted)'])
consumption_fit = consumption["A"]*sym.sin(consumption["w"]*t+consumption["b"])+consumption["c"]

# print(sym.diff(solar_daily["fitfunc"],t))
plt.figure()
plt.title("Sinusodal curve fitting of the electricity consumption")
plt.plot(df['hour'], df['Power (adjusted)'], "ob", label="consumption", linewidth=2)
plt.plot(df['hour'], solar_daily["fitfunc"](tt), "y-", label="fit curve", linewidth=2)
plt.legend(loc="best")
plt.show()
