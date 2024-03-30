clc
clear
close all

t0 = 1; tf = 2; dt = 0.1;
t = t0:dt:tf;

x = -6*(1+t).^(-1) + 3;

plot(t, x, 'LineWidth', 3)
xlabel('time (s)')
ylabel('Optimal x(t)')