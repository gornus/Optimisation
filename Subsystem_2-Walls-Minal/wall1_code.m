
%% initialise
clc
clear all
close all
clear

%% equality constraints
% room dimensions
z=2.2;
y=6;
X= 8;
xmin = 7.2;

% temperature, time and h values
T1= 21;
T2= 15.4;
h1 = 21;
h2 = 10;
%dt = 31540000 % for a year
dt = 1

% shipping container values
outershell_k = 25;
outershell_thickness = 0.0254; 

% total ceiling weight 
solarpanel_mass = X*y*20; %solar pane
ceiling_mass = y*X*6.07; %total mass of just ceiling
mass_load = (ceiling_mass + solarpanel_mass);  %total mass of ceiling + solar panel mass
total_ceiling_force = 9.81*mass_load 
%% ---------OPTIMISING FOR IDEAL VALUES ---------%%

%% initialise 

%set up
x0 = [0.05, 0.05,0.08,0.007,0.05]; %initial point
lb = [0.01, 0.01, 0.01,0.001,0.01]; %adding lower bounds based on manufacturibility and ces data
ub = [];
%adding constraints
A = [0,0,1,1,1];
b = (X-xmin)/2;
Aeq = [];
beq = [];
nonlcon = [];
fun = @(x)(2*z*y)*(T1-T2)*dt/((x(3)/x(1))+(x(4)/x(2))+(x(5)/x(1))+(1/h1)+(1/h2)+(outershell_thickness/outershell_k));

% fmincon ip
options_basic = optimoptions('fmincon','Display','iter','Algorithm','interior-point');
tic; %start timer
IPresult_basic = fmincon(fun,x0,A,b,Aeq,beq,lb,ub, nonlcon, options_basic);
time_fmincon_ip_basic = toc; % end timer

% fmincon SQP optimisation
options2_basic = optimoptions('fmincon','Display','iter', 'Algorithm','sqp');
tic;
sqpresult_basic = fmincon(fun,x0,A,b,Aeq,beq,lb,ub, nonlcon, options2_basic);
time_fmincon_sqp_basic = toc;

% fmincon active set
options3_basic = optimoptions('fmincon','Display','iter','Algorithm','active-set');
tic;
activesetresult_basic = fmincon(fun,x0,A,b,Aeq,beq,lb,ub, nonlcon, options3_basic);
time_fmincon_activeset_basic = toc;

%% -------------OPTIMISING CONSIDERING EACH MATERIAL --------------%%

%% initialise 
%constraints
ub = [];
A = [1,1,1];
b = (X-xmin)/2; 
Aeq = [];
beq = [];
nonlcon = [];

% empty arrays
combo_data = [];
name_data = [];
price_density_data = [];
cost_combo = [];
combo_data_sqp = [];
combo_data_as = [];
combo_data_global = [];
combo_data_ga = [];

%% Reading database of appropriate materials chosen from CES
Mat1_database = readtable('mat_1_database.txt');  
mat1_recordnumber = size(Mat1_database);
Mat2_database = readtable('mat_2_database_final.txt');  
mat2_recordnumber = size(Mat2_database);

m1_T = table2array((Mat1_database(:,2:5)));
m1_conductivity = m1_T(:,1);
m1_cost = m1_T(:,2); %work out for each thickness-- cost per kg
m1_density = m1_T(:,3);
m1_yieldstress = m1_T(:,4) %compressive strength

m2_T = table2array((Mat2_database(:,2:4)));
m2_conductivity = m2_T(:,1);
m2_cost = m2_T(:,2); 
m2_density = m2_T(:,3);

%% collecting other material data for each combination 
for i = 1:mat1_recordnumber(1);
    for j = 1:mat2_recordnumber(1);
        %names
        namem1= table2array(Mat1_database(i,1));
        namem2 = table2array(Mat2_database(j,1));
        names = [namem1 namem2];
        name_data = [name_data; names];        
        %price and density
        m1_cost = m1_T(i,2);
        m2_cost = m2_T(j,2);
        m1_density = m1_T(i,3);
        m2_density = m2_T(j,3);     
        pd_data = [m1_cost m2_cost m1_density m2_density];
        price_density_data = [price_density_data; pd_data];             
   end
end

%% IP and SQP
% finding ideal thicknesses for each combination of materials
for i = 1:mat1_recordnumber(1);
    for j = 1:mat2_recordnumber(1);        
        %% create combinations 
        k_1 = m1_conductivity(i);
        k_2 = m2_conductivity(j);             
        %% min thickness of L1 based on compressive strength
        min_thickness = (total_ceiling_force*0.5)/( 2 * y * m1_yieldstress(i)); %multiplied by 0.5 as wall1 supports 50% of weight 
        lb = [min_thickness,0.001,0.01] ;         
        x0 = [min_thickness+0.05,0.1,0.1]  ;    
        %% function
        fun = @(x)(2*z*y)*(T1-T2)*dt/((x(1)/k_1)+(x(2)/k_2)+(x(3)/k_1)+(1/h1)+(1/h2)+(outershell_thickness/outershell_k));       
        %% running IP and SQP        
        % fmincon interior point optimisation
        options1 = optimoptions('fmincon','Display','iter','Algorithm','interior-point');
        tic; 
        IPresult = fmincon(fun,x0,A,b,Aeq,beq,lb,ub, nonlcon, options1);
        time_fmincon_ip = toc; % end timer             
        x0 = IPresult;  %feed into sqp              
        % fmincon SQP optimisation
        options2 = optimoptions('fmincon','Display','iter', 'Algorithm','sqp');
        tic;
        sqpresult = fmincon(fun,x0,A,b,Aeq,beq,lb,ub, nonlcon, options2);
        time_fmincon_sqp = toc;
        %% append results 
        opt_thick_combo = sqpresult ;
        q_combo = fun(opt_thick_combo);
        new_combo_data = [k_1 k_2 opt_thick_combo q_combo];
        combo_data_sqp = [combo_data_sqp; new_combo_data];                           
   end
end
% finding the optimum combination from results
[q_opt_sqp, I_sqp] = min(combo_data_sqp(:,6));
opt_sqp_combo_data = combo_data_sqp(I_sqp,:); 

%% IP and activeset
% finding ideal thicknesses for each combination of materials
for i = 1:mat1_recordnumber(1);
    for j = 1:mat2_recordnumber(1);       
        %% create combinations 
        k_1 = m1_conductivity(i);
        k_2 = m2_conductivity(j);             
        %% min thickness of L1 based on compressive strength
        min_thickness = (total_ceiling_force*0.5)/( 2 * y * m1_yieldstress(i)); %change for wall2
        lb = [min_thickness,0.001,0.01] ;         
        x0 = [min_thickness+0.05,0.1,0.1]       ;
        %% function
        fun = @(x)(2*z*y)*(T1-T2)*dt/((x(1)/k_1)+(x(2)/k_2)+(x(3)/k_1)+(1/h1)+(1/h2)+(outershell_thickness/outershell_k));        
        %% running IP and active set        
        % fmincon interior point optimisation
        options1 = optimoptions('fmincon','Display','iter','Algorithm','interior-point');
        tic; 
        IPresult = fmincon(fun,x0,A,b,Aeq,beq,lb,ub, nonlcon, options1);
        time_fmincon_ip = toc; % end timer             
        x0 = IPresult;  %feed into active set             
        % fmincon active set
        options3 = optimoptions('fmincon','Display','iter','Algorithm','active-set');
        tic;
        activesetresult = fmincon(fun,x0,A,b,Aeq,beq,lb,ub, nonlcon, options3);
        time_fmincon_activeset = toc;         %end timer
        %% append results 
        opt_thick_combo = activesetresult;
        q_combo = fun(opt_thick_combo);
        new_combo_data = [k_1 k_2 opt_thick_combo q_combo];
        combo_data_as = [combo_data_as; new_combo_data]              
   end
end
%finding optimum combination from results
[q_opt_as, I_as] = min(combo_data_as(:,6));
opt_activeset_combo_data = combo_data_as(I_as,:); 

%% Genetic Algorithm 
% finding ideal thicknesses for each combination of materials
for i = 1:mat1_recordnumber(1);
    for j = 1:mat2_recordnumber(1);       
        %% create combinations 
        k_1 = m1_conductivity(i);
        k_2 = m2_conductivity(j);             
        %% min thickness of L1 based on compressive strength 
        min_thickness = (total_ceiling_force*0.5)/( 2 * y * m1_yieldstress(i)); %change for wall2
        lb = [min_thickness,0.001,0.01] ;         
        x0 = [min_thickness+0.05,0.1,0.1];       
        %% function
        fun = @(x)(2*z*y)*(T1-T2)*dt/((x(1)/k_1)+(x(2)/k_2)+(x(3)/k_1)+(1/h1)+(1/h2)+(outershell_thickness/outershell_k));        
        %% running genetic         
        options6 = optimoptions('ga','ConstraintTolerance',1e-6,'PlotFcn', @gaplotbestf);
        tic; %timer 
        GA_result = ga(fun,3,A,b,Aeq,beq,lb,ub,nonlcon,options6);
        time_ga = toc; %end timer       
        %% append results 
        opt_thick_combo = GA_result ;
        q_combo = fun(opt_thick_combo);
        new_combo_data = [k_1 k_2 opt_thick_combo q_combo];
        combo_data_ga = [combo_data_ga; new_combo_data]   ;                      
   end
end
% find optimum combination
[q_opt_ga, I_ga] = min(combo_data_ga(:,6));
opt_ga_combo_data = combo_data_ga(I_ga,:); 

%% optimum combination results for each solver 
% answer in the format of [k1, k2, l1, l2, l3, heatloss]
opt_ga_combo_data
opt_sqp_combo_data
opt_activeset_combo_data

%% multiobjective  
%calculate cost
combo_data = combo_data_sqp;
row = size(combo_data);
row = row(1);
for v = 1:row
    tot_cost_l1 = z*y*combo_data(v,3)*price_density_data(v,3)*price_density_data(v,1) ;%vol x density x cost
    tot_cost_l2 = z*y*combo_data(v,4)*price_density_data(v,4)*price_density_data(v,2) ;%vol x density x cost
    tot_cost_l3 = z*y*combo_data(v,5)*price_density_data(v,3)*price_density_data(v,1); %vol x density x cost
    totalwall1 = 2*(tot_cost_l1+tot_cost_l2+tot_cost_l3); %both sides of wall
    cost_combo = [cost_combo; totalwall1];
end
% multiobjective graph
figure
x_axis = combo_data(:,6);
y_axis = cost_combo;
index = cellstr(num2str([1:length(cost_combo)]'));
dx = 0.3; dy = 0.3;
scatter(x_axis, y_axis,'filled');
hold on
text(x_axis, y_axis,index);
hold on
edge_bound = boundary(x_axis, y_axis, 0); %creating boundary line with shrink factor of 0
plot(x_axis(edge_bound), y_axis(edge_bound))
grid on
hold on 
% creating a line for y = x as its a 50:50 weighting
thexvalue = linspace(6,16,2);
theyvalue = -thexvalue;
plot(thexvalue,theyvalue)
axis tight
ylabel('Cost (£)')
xlabel('Heat loss (J)') 
zlabel ('Combination number')
title('Multi-Objective Optimisation')

%% Final results and improvement

%choose value based on visual inspection of scatter plot and boundary line
prompt = 'What is the chosen combination number?';
chosencombo = input(prompt) 

% % improvement from ideal values of k and thickness
% q_basic = [fun(IPresult_basic) fun(sqpresult_basic) fun(activesetresult_basic)]; 
% q_basic_average = q_basic (2);  %using sqp result
% q_improvement_from_ideal_perc= ((((combo_data(chosencombo,6)) - (q_basic_average))/(q_basic_average)))*100

% improvement from no insulation
q_no_insulation = ((2*z*y)*(T1-T2))/(outershell_thickness/outershell_k);
q_improvement_from_no_insulation_perc= ((((combo_data(chosencombo,6)) - (q_no_insulation))/(q_no_insulation)))*100

% final selection and results
Material_1 = name_data(chosencombo,1)
Material_2 = name_data(chosencombo,2)
Material_1_and_2_k_values = combo_data(chosencombo,1:2)
L1_L2_L3_values = combo_data(chosencombo,3:5)
wall1_heat_loss = combo_data(chosencombo,6)
wall1_cost = cost_combo(chosencombo)

