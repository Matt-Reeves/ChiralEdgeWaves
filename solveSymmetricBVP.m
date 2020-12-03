clear
%close all
clc

%f = 10; g = f+18.3;
%f = 30; g = f+51.64;

%f = 200; g = f+343.65636615;
f = 300; g = f+515.484548544;
options = odeset('RelTol',1e-6,'AbsTol',1e-6)
[x,y] = ode45(@(t,y) symmetricequilira(t,y,f,g),[0:0.01:10],[1 0],options)


y = y(:,1);
n = exp(y);
dx = x(2)-x(1);
C = -1./(2*pi*sum(x.*y)*dx)
%C = 1./(2*pi*sum(x.*n)*dx)
n = exp(y)*C
%n = n*C 
n = n./(2*pi*trapz(x,x.*n));

M = 2*pi*trapz(x.^3.*n,x)

beta = f/4/pi/C
omega = g/4/beta

figure(10)
plot(x,n*2*pi)