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
import time
from mpl_toolkits.mplot3d import Axes3D

## Data parameters
data_path = "Data\\" # path of the data folder

###############################################################################
###############################################################################
######## Part 1. Formulation of the problem ########

### Model Parameters
# most parameters are here, some are nested in utils as they are not usually changed
resolution = 0.05 # calculation and plotting parameter, smaller = more accurate and takes longer
num_days = 365 # number of days taken account of in this calculation
latitude = np.deg2rad(50.41)
t = sym.Symbol('t') # time parameter for the optimisation problem
length = 6 # length of the 
width = 8
EUR = 0.84 # exchange rate

###############################################################################
### Part 1.1 Setting up the solar irradiance data
solar = Solar_irradiation(latitude, np.arange(1,num_days*24, resolution), resolution)

y_data = solar.get_data()/1000 #converting from Wm^-2 to kWm^-2
plt.figure()
plt.title("Solar Irradiance across the year")
plt.plot((np.arange(1,num_days*24, resolution))/24, y_data, "bo")
plt.xlabel("day")
plt.ylabel("Incident Irradiance $I_{incident} (kWm^{-2})$")
plt.show()

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

###############################################################################
### Part 1.2 Setting up the consumption data 

# Using dataset taken online, with the power consumption adjusted assuming 
# the household in study is average and varies with the national average
# Further adjusted by the region defined using the 2017 consumption data
consmption_df = pd.read_csv(data_path+"energy_consumption.csv")
# converting the time scale from days to hours (to be consistent from above)
consmption_df['hour']=consmption_df['Day']*24

## fit a surrogate model onto the daily average to simplify the algorithm
tt = np.arange(0, num_days*5, 0.1)
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
print("consumption fit: ", consumption_fit)
print("solar daily fit: ", solar_energy_fit)

###############################################################################
### Part 1.3 Importing formulating the solar panel problem

## Reading of the data from the data file
solar_panel_df = pd.read_csv(data_path+"solar_panels.csv")
# Unit conversion to make it consistent to the rest of the formulation
solar_panel_df["Power rating (kW)"] = solar_panel_df["Power rating (W)"]/1000
solar_panel_df["Area (m^2)"] = solar_panel_df["Length (m)"]*solar_panel_df["Width (m)"]
solar_panel_df["Price"] = solar_panel_df["Price (EUR)"]*EUR/30
solar_panel_df["Eff*A"] = solar_panel_df["Module Efficiency"]*solar_panel_df["Area (m^2)"]
# Cleaning up of the dataframe
del solar_panel_df["Power rating (W)"]
del solar_panel_df["OCV (V)"]
del solar_panel_df["SCC (A)"]
del solar_panel_df["Cells"]
del solar_panel_df["Thickness (m)"]
del solar_panel_df["Price (EUR)"]
del solar_panel_df["Annual Degradation (%)"]
del solar_panel_df["NOTC (C)"]
del solar_panel_df["P temp (%/C)"]

# the constraints defined and justified in the report
min_ang = np.deg2rad(20)
max_ang = np.deg2rad(70)
roof_angle = np.arange(min_ang, max_ang, np.deg2rad(1))
# Just to get a better understanding of the problem space and solar panels
print("maximum length of solar panel: ", solar_panel_df["Length (m)"].max())
print("minimum length of solar panel: ", solar_panel_df["Length (m)"].min())
print("maximum roof side length: ", max(roof_width(width, roof_angle)))
print("minimum roof side length: ", min(roof_width(width, roof_angle)))

## Processig of the solar panel data
# any solar panels not in this set are not as cost effective
pareto = pareto_set(solar_panel_df["Name"],solar_panel_df["Price"], \
                    solar_panel_df["Eff*A"],min1=True, min2=False)
# slice this dataframe from the original dataframe
pareto_df = solar_panel_df.loc[solar_panel_df['Name'].isin(pareto[1])]
pareto_df = pareto_df.sort_values(by=["Price"])

# plotting this identified pareto front
plt.figure()
plt.plot(solar_panel_df["Price"],solar_panel_df["Eff*A"],"bo", label = "solar panels")
plt.plot(pareto_df["Price"],pareto_df["Eff*A"],"r-", label = "pareto front", linewidth=1)
plt.plot(pareto_df["Price"],pareto_df["Eff*A"],"yo", label = "pareto set")
plt.legend(loc="best")
plt.title("Power Coefficient vs Equivalent Annual Costs of solar panels")
plt.xlabel("Equivalent Annual Costs (£)")
plt.ylabel("Power Coefficient (kW/panel unit)")
plt.show()



x_data = pareto_df["Price"]
y_data = pareto_df["Eff*A"]
# Data set Needed to be expanded for regression determining the relationship between price and power generation
for i in np.arange(1,4):
    x_data = x_data.append(pareto_df["Price"]*i)
    y_data = y_data.append(pareto_df["Eff*A"]*i)

pareto_front = fit_lin(x_data,y_data)
pareto_fit = lambda x: pareto_front["a"]*x+pareto_front["b"]

# plotting this fit of the pareto front
plt.figure()
plt.plot(x_data, y_data, 'bo')
plt.plot(x_data, pareto_front["fitfunc"](x_data), 'r-')
plt.title("Relationship between Price and Power Generation")
plt.xlabel("Equivalent Annual Costs (£)")
plt.ylabel("Power Coefficient (kW)")
plt.show()

###############################################################################
###############################################################################
######## Part 2. Exploration of the problem space ########

### Part 2.1 Explore how the change in roof angle affect the objective function

consumption_fit = lambda t: 1.30338864562213 * np.sin(0.0172004453662375*t + 1.62737939447036) + 11.9623058930769
solar_fit = lambda t: 0.2*(5.79292108750245 * np.sin(0.0162900269231784*t - 1.2108983113935) + 9.92337227046363)

beta = np.arange(min_ang, max_ang, np.deg2rad(2))
gen_rate = 3
panel_costs = []
energy_costs = []
total_costs = []


for i in beta:
    generation = lambda t: gen_rate*np.sin(np.pi/2-np.deg2rad(latitude)+solar.delta(t)+i)*solar_fit(t)
    energy_cost = bills(generation,consumption["fitfunc"],[0, 365])
    panel_cost = pareto_fit(gen_rate)
    energy_costs.append(energy_cost)
    panel_costs.append(panel_cost)
    total_costs.append(energy_cost+panel_cost)

plt.figure()
plt.plot(np.rad2deg(beta), energy_costs, 'r-', label = "total costs")
#plt.plot(np.rad2deg(beta), total_costs, 'b-')
min_index = total_costs.index(min(total_costs))
plt.plot(np.rad2deg(beta[min_index]), min(total_costs), 'bo')
plt.xlabel("Roof angle beta (rads)")
plt.ylabel("equivalent annual costs")
plt.title("study on the effect of roof angle")
print("minimum point beta: ", np.rad2deg(beta[min_index]))

### Part 2.2 Explore how the change in generation rate affect the objective function

beta = np.deg2rad(20)
gen_rate = np.arange(0.24, 10, 0.05)
panel_costs = []
energy_costs = []
total_costs = []

for i in gen_rate:
    generation = lambda t: i*np.sin(np.pi/2-np.deg2rad(latitude)+solar.delta(t)+beta)*solar_fit(t)
    energy_cost = bills(generation,consumption["fitfunc"],[0, 365])
    panel_cost = pareto_fit(i)
    energy_costs.append(energy_cost)
    panel_costs.append(panel_cost)
    total_costs.append(energy_cost+panel_cost)

plt.figure()
plt.plot(gen_rate, total_costs, 'r-', label = "total costs")
#plt.plot(np.rad2deg(beta), total_costs, 'b-')
min_index = total_costs.index(min(total_costs))
plt.plot(gen_rate[min_index], min(total_costs), 'bo')
plt.xlabel("Power Generation (kW)")
plt.ylabel("equivalent annual costs")
plt.title("study on the effect of Power Generation rate")
print("minimum point gen rate: ", gen_rate[min_index])

###############################################################################
###############################################################################
######## Part 3. Optimisation ########

#### Part 3.1 Using the fit function and contour ####
beta = np.arange(np.deg2rad(20), np.deg2rad(30), np.deg2rad(0.1))
gen_rate = np.arange(3, 6, 0.02)
roof_angs = np.array([])
gen_rates = np.array([])
total_costs = np.zeros((len(beta),len(gen_rate)))

for i in range(len(beta)):
    for j in range(len(gen_rate)):
        generation = lambda t: gen_rate[j]*np.sin(np.pi/2-np.deg2rad(latitude)+solar.delta(t)+beta[i])*solar_fit(t)
        energy_cost = bills(generation,consumption_fit,[0, 365])
        panel_cost = pareto_fit(gen_rate[j])#+roof_width(width, beta[i])*length*roof_unit_costs/roof_life
        gen_rates = np.append(gen_rates,gen_rate[j])
        roof_angs = np.append(roof_angs,beta[i])
        total_costs[i,j]=energy_cost+panel_cost


#plt.figure()
#plt.contour(gen_rates, roof_angs, total_costs)
#
##ax = fig.add_subplot(111, projection='3d')
##
##ax.scatter(gen_rates, roof_angs, total_costs, c='r', marker='o')
#
#	
## Find index of minimum value from 2D numpy array
#result = np.where(total_costs == np.amin(total_costs))
#plt.plot(gen_rates[result[0]], roof_angs[result[1]], np.amin(total_costs), 'bo')
#
#plt.xlabel('Generation Rates ($kWm^-2$)')
#plt.ylabel('Roof angle (degree)')
#plt.title("Problem Space Contour")
#plt.show()

print("minimum cost: %.2f; angle: %.2f; gen rate: %.2f" %(np.amin(total_costs), beta[result[1]], gen_rate[result[0]]))

# 
target = 3
for i in range(pareto_df.shape[0]):
    number = np.ceil(target/pareto_df.iloc[i]["Eff*A"])
    cost = number*pareto_df.iloc[i]["Price"]
    print("model: %s; cost: %.2f; number: %d" %(pareto_df.iloc[i]["Name"],cost, number))
    
print("\n\n")

#### Part 3.2 Direct Search ####
overall_start=time.time()
beta = np.arange(np.deg2rad(20), np.deg2rad(25), np.deg2rad(0.1))

print("Starting Direct Search")

for i in range(pareto_df.shape[0]):
    no_panels = np.array([])
    roof_angs = np.array([])
    total_costs = np.array([])

    start_time = time.time()
    for j in beta:
        # loop through the beta
        max_panels = max_tiling(pareto_df.iloc[i]["Length (m)"],pareto_df.iloc[i]["Width (m)"],\
                               length, roof_width(width, j))
#         print("Model: %s; roof angle: %.2f; max panels: %d" %(pareto_df.iloc[i]["Name"],np.rad2deg(j),max_panels))

        for k in range(int(max_panels)):
            generation = lambda t: k*pareto_df.iloc[i]["Eff*A"]*np.sin(np.pi/2-np.deg2rad(latitude)+solar.delta(t)+j)*solar_fit(t)
            energy_cost = bills(generation,consumption_fit,[0, 365])
            panel_cost = k*pareto_df.iloc[i]["Price"]#+roof_width(width, j)*length*roof_unit_costs/roof_life
            no_panels = np.append(no_panels,k)
            roof_angs = np.append(roof_angs,j)
            total_costs = np.append(total_costs,energy_cost+panel_cost)
#             print("No. Panels: %d; energy cost: %.2f; panel cost: %.2f" %(k,energy_cost,panel_cost))


    index = total_costs.argsort()[::1]
    total_costs = total_costs[index] #sort from largest to smallest
    no_panels = no_panels[index]
    roof_angs = roof_angs[index]
    print("Model: %s; minimum cost: %.2f; no. of panels: %d; roof angle: %.2f" 
          %(pareto_df.iloc[i]["Name"],total_costs[0], no_panels[0],np.rad2deg(roof_angs[0])))
    
    print("--- loop %d: %s seconds ---" % (i, time.time() - start_time))
print("--- total: %s seconds ---" % (time.time() - overall_start))    