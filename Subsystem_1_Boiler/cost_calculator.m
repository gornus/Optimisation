function [cost] = cost_calculator(x, i, j)
  c1 = [2.39,0.27,0.638,17.9,5.09,0.0380];     %cost per kg (£/kg)
  c2 = [3.56,1.75,1.93,4.33,1.68,4.68];        %cost per kg (£/kg)
  
  %cost of material 1
  r2 = 0.223+x(1);         %radius of internal and material 1
  L2 = 0.446+(2*x(1));     %length of internal and material 1
  Vb = pi*r2^2*L2;         %volume of internal and material 1
  Vm1 = Vb-0.07;           %volume of material 1
  Mm1 = x(3)*Vm1;          %mass of material 1
  cost1 = c1(i)*Mm1;
  
  %cost of material 2
  r3 = r2+x(2);            %radius of entire boiler
  L3 = L2+(2*x(2));        %length of entire boiler
  Ve = pi*r3^2*L3;         %volume of boiler
  Vm2 = Ve-Vb;             %volume of material 2
  Mm2 = x(4)*Vm2;          %mass of material 2
  cost2 = c2(j)*Mm2;
  
  cost = cost1+cost2;
end