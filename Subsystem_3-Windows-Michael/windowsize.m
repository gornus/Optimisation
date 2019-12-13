DeltaT = 5.6;
h0 = 10;
hi = 30;
kg = 1;
ka = 0.016;
lg = 0.01;
la = 0.01;
A = 0:.1:13.2;
Q=(((DeltaT))*A)/((1/h0)+(1/hi)+(2*lg/kg)+(la/ka));
figure
plot(A,Q);
p=polyfit(A,Q,1)
title('Heating Energy');
ylabel('Energy Loss, W') 
xlabel('Area of window,m^2')