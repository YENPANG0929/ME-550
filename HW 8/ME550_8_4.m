clc
clear
close all

syms c1 c2
f1 = c1 + c2;
f2 = c1*exp(1) + c2*exp(-1);

[c1, c2] = solve([f1 == 1, f2 == 1], [c1, c2]);
C1 = double(c1)
C2 = double(c2)

t0 = 0; tf = 1; dt = 0.01;
t = t0:dt:tf;

x = c1*exp(t) + c2*exp(-t);
x_dot = c1*exp(t) - c2*exp(-t);
u = x_dot + x;
figure(1)
plot(t, x, 'LineWidth', 3)
xlabel('time (s)')
ylabel('Optimal x(t)')
figure(2)
plot(t, u, 'LineWidth', 3)
xlabel('time (s)')
ylabel('Optimal u(t)')
