function [c,ceq] = nlcon(x)
  %mass of material 1
  r2 = 0.2673+x(1);         %radius of internal and material 1
  L2 = 0.5346+(2*x(1));     %length of internal and material 1
  Vb = pi*r2^2*L2;         %volume of internal and material 1
  Vm1 = Vb-0.12;           %volume of material 1
  Mm1 = x(3)*Vm1;
  
  %mass of material 2
  r3 = r2+x(2);            %radius of entire boiler
  L3 = L2+(2*x(2));        %length of entire boiler
  Ve = pi*r3^2*L3;         %volume of boiler
  Vm2 = Ve-Vb;             %volume of material 2
  Mm2 = x(4)*Vm2;
  
  %combined mass
  Mb = Mm1+Mm2+70;   %total mass of boiler
  c = Mb - 250;      %maximum mass of boiler
  ceq = [];          %empty constraints
end