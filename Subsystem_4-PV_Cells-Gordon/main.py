############################
# Created on Wed Dec  4 15:51:29 2019
# author: Gordon Cheung Yat Hei
# CID: 01083012
# Project: Optimisation
############################

import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
from solar_irradiation import Solar_irradiation
from utils import *
import pandas as pd

## Data parameters
data_path = "Data\\" # path of the data folder

######## Part 1. Formulation of the problem ########

### Model Parameters
resolution = 0.05
num_days = 365 # number of days taken account of in this calculation
latitude = np.deg2rad(50.41)
t = sym.Symbol('t') # time parameter for the optimisation problem

### Part 1.1 Setting up the solar irradiance data
solar = Solar_irradiation(latitude, np.arange(1,num_days*24, resolution), resolution)

y_data = solar.get_data()/1000 #converting from Wm^-2 to kWm^-2
#plt.figure()
#plt.title("Solar Irradiance across the year")
#plt.plot((np.arange(1,num_days*24, resolution))/24, y_data, "bo")
#plt.xlabel("day")
#plt.ylabel("Incident Irradiance $I_{incident} (kWm^{-2})$")
#plt.show()

## get the daily averages of solar irradiation
yavg = solar.get_daily_average()/1000 #converting from Wm^-2 to kWm^-2


## fit a surrogate model onto the daily average to simplify the algorithm
day = np.arange(num_days)
solar_daily = fit_sin(day, yavg) #fitting based in days
# symbolic function of the fit
solar_daily_fit = solar_daily["A"]*sym.sin(solar_daily["w"]*(t)+solar_daily["b"])+solar_daily["c"]
solar_energy_fit = solar_daily_fit*24
# unit:  time = day; irradiance = kWm^-2
plt.figure()
plt.title("Incident Irradiance")
plt.plot((np.arange(1,num_days*24, resolution))/24, y_data, "bo", label = "incident irradiance")
plt.plot(day, yavg, "ro", label="daily average incident irradiance", linewidth=0.5)
plt.plot(day, solar_daily["fitfunc"](day), "y-", label="sin fit for the daily average", linewidth=2) # based in hours
plt.legend(loc=8)
plt.xlabel("day")
plt.ylabel("Incident Irradiance $I_{incident} (kWm^{-2})$")
plt.show()



## compute the maximum and minimum days for later optimisation
max_day = np.argmax(yavg)
min_day = np.argmin(yavg)
yId = solar.get_data()/1000
print("max gen day: ", max_day)
print("min gen day: ", min_day)
# Using the maximum and minimum to slice the days out for later optimisation
max_data = yId[int(max_day*24/resolution):int((max_day*24+23)/resolution)]
min_data = yId[int(min_day*24/resolution):int((min_day*24+23)/resolution)]
#min_data = yId[int(min_day/resolution):int((min_day+23)/resolution)]
plt.figure()
plt.title("Maximum and Minimum Incident Irradiance Days")
plt.plot(np.arange(0,23,resolution),max_data, "y-", label="max day", linewidth=1)
plt.plot(np.arange(0,23,resolution),min_data, "r-", label="min day", linewidth=1)
plt.xlabel("hour")
plt.ylabel("Incident Irradiance $I_{incident} (kWm^{-2})$")
plt.legend(loc="best")
plt.show()


## fit a surrogate model onto the daily average to simplify the algorithm
tt = np.arange(0, 23,resolution)
solar_max = fit_sin(tt, max_data)
solar_min = fit_sin(tt, min_data)


### Part 1.2 Setting up the consumption data 

# Using dataset taken online, with the power consumption adjusted assuming 
# the household in study is average and varies with the national average
# Further adjusted by the region defined using the 2017 consumption data
consmption_df = pd.read_csv(data_path+"energy_consumption.csv")
# converting the time scale from days to hours (to be consistent from above)
consmption_df['hour']=consmption_df['Day']*24

## fit a surrogate model onto the daily average to simplify the algorithm
tt = np.arange(0, num_days, 0.1)
consumption = fit_sin(consmption_df['Day'], consmption_df['Power (adjusted)']/(30)) # scaled to daily average instead of monthly total
# symbolic function of the fit
consumption_fit = consumption["A"]*sym.sin(consumption["w"]*t+consumption["b"])+consumption["c"]
# units:  time = day; power = kWh

## Plotting the curve of daily average solar irradiation with daily average consumption
plt.figure()
plt.title("Estimated electricity consumption per day")
plt.plot(consmption_df['Day'], consmption_df['Power (adjusted)']/(30), "yo", label="daily consumption", linewidth=2) # based in hours
plt.plot(tt, consumption["fitfunc"](tt), "b-", label="daily consumption fit curve", linewidth=1)
plt.legend(loc="best")
plt.xlabel("day")
plt.ylabel("Daily Energy Consumption $E_{C} (kWh)$")
plt.show()
#

max_day = np.argmax(consumption["fitfunc"](tt))
min_day = np.argmin(consumption["fitfunc"](tt))
print("max use day: ", max_day)
print("min use day: ", min_day)


## verifying the converted function matches the data
year = (365)
total = sym.integrate(consumption_fit, (t, 0, year))
print (total)
print(consumption_fit)
print(solar_energy_fit)

### Part 1.3 Importing formulating the solar panel problem



