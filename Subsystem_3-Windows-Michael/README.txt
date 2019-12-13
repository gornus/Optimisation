READ ME
Main: Contains the main FMINCON optimisation within MATLAB which combines the simulated light equations with the central heat loss equation. Outputting the optimal size of the window, balancing lighting and heating energy requirements. Obj, Nonlcon, feed into the Main file for the optimisation. 

Simulation: takes the data from the RAW Simulation data set as a LightEnergy csv and plots this. From this a meta objective function can be found describing the relationship between the window size and the lighting energy.

Window Size & Window Gap: Plots of the objective function and changing various variables. The upper and lower bounds are set within these functions. For Window Gap, this can also be used for the Window Materials for analytical analysis of the various factors.




