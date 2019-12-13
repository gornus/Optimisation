close all
clear all

%% Material Properties

%Inner Material Properties
%Stainless,Iron,Alloy-Steel,Titanium,Bronze,Concrete
k1 = [19.5,36.5,43.6,7.50,50,1.4];           %thermal conductivity value (W/mK)
d1 = [7740,7100,7800,4610,8375,2400];        %density (kg/m^3)
c1 = [2.39,0.27,0.638,17.9,5.09,0.0380];     %cost per kg (£/kg)
Oy = [699,438,1030,896,320,1.13];            %yield strength (MPa)

%Insulating Material Properties
%Polyurethane,Butyl,FPF,Cellulose,PS,Cork
k2 = [0.24,0.105,0.042,0.14,0.13,0.044];     %thermal conductivity value (W/mK)
d2 = [1180,930,15.5,990,1.04,197.5];         %density (kg/m^3)
c2 = [3.56,1.75,1.93,4.33,1.68,4.68];        %cost per kg (£/kg)



%% Big Boi Code
%fmincon(fun,x0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
%fun = the function to minimise
%x0 = initial guess
%A,B = constraints
%Aeq,Beq = equality constraints
%LB,UB = lower/ upper bounds
%NONLCON = non-linear constraints

heat_loss = [];
cost = [];
n = 0;

for i = 1:6
    for j = 1:6
        %input material properties
        d10 = d1(i); d20 = d2(j);
        k10 = k1(i); k20 = k2(j);
        OyO = Oy(i);
                
        %load things for fmincon 
        fx = @(x) objective(x);
        x0 = [0.01,0.3,d10,d20,k10,k20];  %initial guess
        A = [1,1,0,0,0,0;   %linear constraint variables
             1,0,0,0,0,0];
        B = [0.277;         %linear constraint values
             1.68/OyO];
        Aeq = [0,0,1,0,0,0;    %set all material properties
               0,0,0,1,0,0;
               0,0,0,0,1,0;
               0,0,0,0,0,1];
        Beq = [d10;            %material properties from table
               d20;
               k10;
               k20];
        lb = [0.001,0.001,0,0,0,0]; ub = [];  %upper and lower bounds
        nonlcon = @nlcon;
        
        %fmincon time bois
        y = fmincon(fx,x0,A,B,Aeq,Beq,lb,ub,nonlcon);
        
        %store values for plotting
        n = n+1;
        heat_loss(n) = objective(y);          %heat loss
        cost(n) = cost_calculator(y,i,j);     %cost of materials
        
        
    end
end



%% Plot Results

figure      %all material combination plots
c = [0,0,0,0,0,0,2,2,2,2,2,2,4,4,4,4,4,4,6,6,6,6,6,6,8,8,8,8,8,8,10,10,10,10,10,10]; %colours for plot
hold on
scatter(heat_loss,cost,80,c,'filled')  %plot data
scatter(157,454,80,0,'filled','d')    %plot benchmark
hold off
title('All Material Combinations Comparison')
xlabel('Heat Loss (J/hr)')
ylabel('Cost (£)')
legend('optimised solutions','benchmark')
grid on

%filtered material combination plots
heat_loss_filt = [];
cost_filt = [];
for i = 1:6
    j = (i*6) - 3;   %third insulating material (FPF)
    heat_loss_filt(i) = heat_loss(j);
    cost_filt(i) = cost(j);
end
c_filt = [0,2,4,6,8,10]; %colours for filtered plot
figure
scatter(heat_loss_filt,cost_filt,80,c_filt,'filled')
title('FPF Insulator Material Combinations Comparison')
xlabel('Heat Loss (J/hr)')
ylabel('Cost (£)')
grid on


