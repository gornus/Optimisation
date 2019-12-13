function [fx] = objective(x)
  %Variables
  T0 = 60; T4 = 21;   %temperatures
  hw = 700; ha = 21;  %convection coefficents
  r1 = 0.2673;         %inner radius (m) - solved analytically
  L1 = 0.5346;         %inner length (m) - solved analytically
  
  %Equations
  r2 = r1+x(1);          %middle radius
  r3 = r2+x(2);          %outer radius
  L2 = L1+x(1);          %middle length
  L3 = L2+x(2);          %outer length
  Ai = 2*pi*r1*(r1+L1);  %internal
  Ae = 2*pi*r3*(r3+L3);  %external

  %Convection
  Ei = 1/(hw*Ai);  %internal
  Ee = 1/(ha*Ae);  %external

  %Conduction - k1=x(1) , k2=x(2)
  Ec1 = (x(1)*log(r2/r1)) / (2*pi*x(5)*((L1*x(1))+(r1^2*log(r2/r1))));  %structure
  Ec2 = (x(2)*log(r3/r2)) / (2*pi*x(6)*((L2*x(2))+(r2^2*log(r3/r2))));  %insulator

  fx = (T0-T4)/(Ei+Ec1+Ec2+Ee); % minimise
end