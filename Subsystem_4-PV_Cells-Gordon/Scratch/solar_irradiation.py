############################
# Created on Wed Dec  4 15:50:24 2019
# author: Gordon Cheung Yat Hei
# CID: 01083012
# Project: Optimisation
############################

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize

class Solar_irradiation:
    # Class encompassing all the functions used to calculate the solar irradiation data 
    # from first principles
    def __init__(self, phi, t, resolution):
        self.phi = phi #local latititude
        self.t = t #time array, in hours
#        self.yId = yId
#        self.yavg = yavg
        self.resolution = resolution #resolution in hours
        
        ## computing the data arrays of local irradiation and daily average
        # initialisation of the data arrays
        ti = len(self.t)
        self.yId = np.zeros(ti)
        self.yavg = np.zeros(int(np.floor(ti/(24/self.resolution)))+1) # the array for daily average
        for i in range(ti):
            # calculating the direct solar irradiance
            data = self.Id(self.t[i],self.phi)
            self.yId[i] = data
            index = int(np.floor(i/(24/self.resolution)))-1
            # using the local solar irradiance data to calculate the daily average
            self.yavg[index] = self.yavg[index] + data*self.resolution
        self.yavg = self.yavg / (24)
        
    def delta(self, t):
        return np.deg2rad(23.45)*np.sin(np.deg2rad((360/365)*(t/24-81)))

    # hour angle at the local solar time
    def h(self, t):
        return -np.pi + np.deg2rad(360*t/24)

    # cos of solar zenith angle
    def costheta(self, t, phi):
        return np.sin(phi)*np.sin(self.delta(t))+np.cos(phi)*np.cos(self.delta(t))*np.cos(self.h(t))

    # air mass calculating solar solar irradiance accounting for the atmosphere
    def AM(self, t,phi):
        r = 6371*1000 # in meters
        return np.sqrt((r*self.costheta(t,phi))**2+2*r+1)-r*self.costheta(t,phi)

    # outputting the incident solar irradiance
    def Id(self, t, phi):
        IS = 1353 # in watts per meter^2
        factor = 1.25 # the factor for taking account DHI as a ratio of DNI
        return factor*1.1*IS*0.7**(self.AM(t,phi)**0.678)
    
    # 
    def get_data(self):
        return self.yId
    
    def get_daily_average(self):
        return self.yavg