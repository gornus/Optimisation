data=csvread('lightenergy.csv');
x=data(:,1);
y=data(:,2);
%plot(x,y,'o');
p = polyfit(x,y,3)
x2 = 0:0.1:13.2;
y2 = polyval(p,x2);
figure
plot(x,y,'o',x2,y2)
title('Light Energy');
xlabel('Window Area, m^2');
ylabel('Artifical Lighting Energy, W');