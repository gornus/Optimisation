DeltaT = 5.6;
h0 = 10;
hi = 30;
kg = 1;
lg = 0.01;
ka = 0.016;
A = 10.32;
loss = [];
width = linspace(0,0.03,30);
for la = 0.0001:.001:0.03;
Q=((DeltaT)*A)/((1/h0)+(1/hi)+((2*lg)/kg)+(la/ka)); 
loss = [loss; [Q]];
end
figure
plot(width,loss)
ylabel('Energy Loss, W')
xlabel('size of air gap, m')