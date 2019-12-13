############################
# Created on Thu Dec 12 06:31:40 2019
# author: Gordon Cheung
# CID: 01083012
# Project: Optimisation
############################
import pandas as pd
import sympy as sym
import numpy as np
import math
import matplotlib.pyplot as plt
import time
from mpl_toolkits.mplot3d import Axes3D
from scipy import optimize
from solar_irradiation import Solar_irradiation

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


length = 6
width = 8
EUR = 0.84 # exchange rate

## Reading of the data
solar_panel_df = pd.read_csv("Optimisation\\Subsystem_4-PV_Cells-Gordon\\Data\\solar_panels.csv")
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
#solar_panel_df.head()

def roof_width(width, beta):
    return (width*np.sin(beta))/np.sin(np.pi-2*beta)

min_ang = np.deg2rad(20)
max_ang = np.deg2rad(70)
roof_angle = np.arange(min_ang, max_ang, np.deg2rad(1))

#plt.plot(roof_angle, length*roof_width(width, roof_angle), label = "roof area")
#plt.plot(roof_angle, roof_width(width, roof_angle), label="roof side")
#plt.legend(loc="best")

print("maximum length of solar panel: ", solar_panel_df["Length (m)"].max())
print("minimum length of solar panel: ", solar_panel_df["Length (m)"].min())
print("maximum roof side length: ", max(roof_width(width, roof_angle)))
print("minimum roof side length: ", min(roof_width(width, roof_angle)))

def max_tiling(p_length, p_width, r_length, r_width):
    if r_width >= p_length:
        # if the panels can be fitted perpendicularly
        number = np.floor(r_width/p_length)*np.floor(r_length/p_width)
        if r_width-p_length*np.floor(r_width/p_length) >= p_width:
            # if more panels can be fitted parallely in the remaining space
            number = number + np.floor(r_width-p_length*np.floor(r_width/p_length)/p_width)
    elif r_width >= p_width:
        # otherwise if we can fit the panels parallely
        number = np.floor(r_width/p_width)*np.floor(r_length/p_length)
    else: number = 0
    return number


# do not need symbolic integration as they are daily averages

def bills(energy, consumption, domain):    
    buy_price = 12.5/100 #grid power selling price/kWh
#    sell_price = 5.5/100 #grid power buying price/kWh
    sell_price = 12.5/100 #grid power buying price/kWh
    cost = 0 # initialising costs
    sell = 0
    needed = 0
    for i in np.arange(domain[0],domain[1]):
        difference = consumption(i)-energy(i)
        if difference > 0:
            # if there is more consumption than energy generated
            cost = cost + difference*buy_price
        else:
            # if the excess energy can be sold back to the grid
            sell = sell + difference*sell_price
    return cost-sell

## determine the maximum possible costs of purchasing solar panels
## Max roof area happens when max roof angle
#
#max_roof_width = roof_width(width, max_ang)
#min_roof_width = roof_width(width, min_ang)
#print("max area: ", max_roof_width*length)
#print("min area: ", min_roof_width*length)
#num_models = solar_panel_df.shape[0] #number of solar panels in the data set
#max_tile = []
#min_tile = []
#for i in range(num_models):
##     print("length: ", solar_panel_df.iloc[i]["Length (m)"])
##     print("width: ", solar_panel_df.iloc[i]["Width (m)"])
#    max_tile.append(max_tiling(solar_panel_df.iloc[i]["Length (m)"], \
#                               solar_panel_df.iloc[i]["Width (m)"], length, max_roof_width))
#    min_tile.append(max_tiling(solar_panel_df.iloc[i]["Length (m)"], \
#                               solar_panel_df.iloc[i]["Width (m)"], length, min_roof_width))
## print(max_tile)
#solar_panel_df.insert(solar_panel_df.shape[1],"max tiles", max_tile) # figure out the 
#solar_panel_df.insert(solar_panel_df.shape[1],"min tiles", min_tile) # figure out the 
#solar_panel_df.head()

#plt.figure()
#plt.plot(solar_panel_df["Price"],solar_panel_df["Eff*A"],"bo")

def pareto_set(name, data_1, data_2, min1=True, min2=True):
    # sort by data_1 to find the dominating points
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
#     print(data_1)
#     print("first data: ", data_1[index[0]], " ", data_2[index[0]], " ", name[index[0]])
#     print(index)
#     print(index[len(index)-1])
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

# any solar panels not in this set are not as cost effective
pareto = pareto_set(solar_panel_df["Name"],solar_panel_df["Price"], \
                    solar_panel_df["Eff*A"],min1=True, min2=False)
# print(pareto[0])
# print(pareto[1])

# As there are only three solar panels in the pareto set
# Needed to be expanded for regression determining the relationship between price and power generation
# for panel_model in pareto:
#     print(solar_panel_df[solar_panel_df["Name"]==panel_model])
pareto_df = solar_panel_df.loc[solar_panel_df['Name'].isin(pareto[1])]
pareto_df = pareto_df.sort_values(by=["Price"])
#plt.figure()
#plt.plot(solar_panel_df["Price"],solar_panel_df["Eff*A"],"bo", label = "solar panels")
#plt.plot(pareto_df["Price"],pareto_df["Eff*A"],"r-", label = "pareto front", linewidth=1)
#plt.plot(pareto_df["Price"],pareto_df["Eff*A"],"yo", label = "pareto set")
#plt.legend(loc="best")
#plt.title("Power Coefficient vs Equivalent Annual Costs of solar panels")
#plt.xlabel("Equivalent Annual Costs (£)")
#plt.ylabel("Power Coefficient (kW/panel unit)")
#plt.show()

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

x_data = pareto_df["Price"]
y_data = pareto_df["Eff*A"]
for i in np.arange(1,4):
    x_data = x_data.append(pareto_df["Price"]*i)
    y_data = y_data.append(pareto_df["Eff*A"]*i)
print(x_data.shape)
# print(x_data)
# print(y_data)
pareto_front = fit_lin(x_data,y_data)
# pareto_fit = lambda x: pareto_front["a"]*np.exp(pareto_front["b"]*x)+pareto_front["c"]
pareto_fit = lambda x: pareto_front["a"]*x+pareto_front["b"]
print(pareto_front["a"])
print(pareto_front["b"])
print(pareto_fit)
#plt.figure()
#plt.plot(x_data, y_data, 'bo')
#plt.plot(x_data, pareto_front["fitfunc"](x_data), 'r-')
#plt.title("Relationship between Price and Power Generation")
#plt.xlabel("Equivalent Annual Costs (£)")
#plt.ylabel("Power Coefficient (kW)")
#plt.show()
# not a good fit but worth trying

# consumption_fit = 1.30338864562213*sym.sin(0.0172004453662375*t + 1.62737939447036) + 11.9623058930769
# solar_fit = 5.62575334261198*sin(0.01716050340221*t - 1.34948474433562) + 10.2277234035447


# pareto_fit = 0.0560335267904956*exp(19.7798406889663*x) + 96.2783894507582
# roof_factor = sin(np.pi/2-np.deg2rad(latitude)+delta(t)+beta)
# delta = np.deg2rad(23.45)*np.sin(np.deg2rad((360/365)*(t/24-81)))
consumption_fit = lambda t: 1.30338864562213 * np.sin(0.0172004453662375*t + 1.62737939447036) + 11.9623058930769
solar_fit = lambda t: 0.2*(5.62575334261198 * np.sin(0.01716050340221*t - 1.34948474433562) + 10.2277234035447)
#delta = lambda t: np.deg2rad(23.45)*np.sin(np.deg2rad((360/365)*(t/24-81)))
roof_factor = lambda beta, t: np.sin(np.pi/2-np.deg2rad(latitude)+delta(t)+beta)
#
#    
#beta = np.arange(min_ang, max_ang, np.deg2rad(2))
#latitude = np.deg2rad(50.41)
#gen_rate = 3
#panel_costs = []
#energy_costs = []
#total_costs = []
#for i in beta:
#    generation = lambda t: gen_rate*np.sin(np.pi/2-np.deg2rad(latitude)+solar.delta(t)+i)*solar_fit(t)
#    energy_cost = bills(generation,consumption_fit,[0, 365])
#    panel_cost = pareto_fit(gen_rate)#+roof_width(width, i)*length*roof_unit_costs/roof_life
#    energy_costs.append(energy_cost)
#    panel_costs.append(panel_cost)
#    total_costs.append(energy_cost+panel_cost)
#
#plt.figure()
#plt.plot(np.rad2deg(beta), energy_costs, 'r-', label = "total costs")
##plt.plot(np.rad2deg(beta), total_costs, 'b-')
#min_index = total_costs.index(min(total_costs))
#plt.plot(np.rad2deg(beta[min_index]), min(total_costs), 'bo')
#plt.xlabel("Roof angle beta (rads)")
#plt.ylabel("equivalent annual costs")
#plt.title("study on the effect of roof angle")
#print("minimum point beta: ", np.rad2deg(beta[min_index]))
#
#
#beta = np.deg2rad(20)
#latitude = np.deg2rad(50.41)
#gen_rate = np.arange(0.24, 10, 0.05)
#panel_costs = []
#energy_costs = []
#total_costs = []
#roof_life = 30
#roof_unit_costs = 90
#
#
#for i in gen_rate:
#    generation = lambda t: i*np.sin(np.pi/2-np.deg2rad(latitude)+solar.delta(t)+beta)*solar_fit(t)
#    energy_cost = bills(generation,consumption_fit,[0, 365])
#    panel_cost = pareto_fit(i)+roof_width(width, beta)*length*roof_unit_costs/roof_life
#    energy_costs.append(energy_cost)
#    panel_costs.append(panel_cost)
#    total_costs.append(energy_cost+panel_cost)
#
#plt.figure()
#plt.plot(gen_rate, total_costs, 'r-', label = "total costs")
##plt.plot(np.rad2deg(beta), total_costs, 'b-')
#min_index = total_costs.index(min(total_costs))
#plt.plot(gen_rate[min_index], min(total_costs), 'bo')
#plt.xlabel("Power Generation (kW)")
#plt.ylabel("equivalent annual costs")
#plt.title("study on the effect of Power Generation rate")
#print("minimum point gen rate: ", gen_rate[min_index])


beta = np.arange(np.deg2rad(20), np.deg2rad(30), np.deg2rad(0.1))
gen_rate = np.arange(3, 6, 0.02)
roof_angs = np.array([])
gen_rates = np.array([])
total_costs = np.zeros((len(beta),len(gen_rate)))

for i in range(len(beta)):
    for j in range(len(gen_rate)):
        generation = lambda t: gen_rate[j]*np.sin(np.pi/2-np.deg2rad(latitude)+solar.delta(t)+beta[i])*solar_fit(t)
        energy_cost = bills(generation,consumption_fit,[0, 365])
        panel_cost = pareto_fit(gen_rate[j])+roof_width(width, beta[i])*length*roof_unit_costs/roof_life
        energy_costs.append(energy_cost)
        panel_costs.append(panel_cost)
        total_costs = np.append(total_costs, energy_cost+panel_cost)
        gen_rates = np.append(gen_rates,j)
        roof_angs = np.append(roof_angs,i)
        total_costs[i,j]=energy_cost+panel_cost


plt.figure()
plt.contour(gen_rate, beta, total_costs)
#ax = fig.add_subplot(111, projection='3d')
#
#ax.scatter(gen_rates, roof_angs, total_costs, c='r', marker='o')

	
# Find index of minimum value from 2D numpy array
#result = np.where(total_costs == np.amin(total_costs))
#plt.plot(gen_rate[result[0]], beta[result[1]], np.amin(total_costs), 'bo')


plt.xlabel('Generation Rates ($kWm^-2$)')
plt.ylabel('Roof angle (degree)')
plt.title("Problem Space Contour")
plt.show()

print("minimum cost: %.2f; angle: %.2f; gen rate: %.2f" %(np.amin(total_costs), beta[result[1]], gen_rate[result[0]]))

#plt.plot(gen_rate, energy_costs, 'r-', label = "energy costs")
##plt.plot(np.rad2deg(beta), total_costs, 'b-')
#min_index = energy_costs.index(min(energy_costs))
#plt.plot(gen_rate[min_index], min(energy_costs), 'bo')
#plt.xlabel("Power Generation (kW)")
#plt.ylabel("equivalent annual costs")
#plt.title("study on the effect of Power Generation rate")
#print("minimum point gen rate: ", gen_rate[min_index])