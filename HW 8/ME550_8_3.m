clc
clear
close all

syms c5 c6 c7 c8
f1 = c5 + c6 + c7;
f2 = 16*c5 + 8*c6 + c7;
f3 = c5/5 + c6/4 + c7 + c8;
f4 = 32*c5/5 + 4*c6 + 2*c7 + c8;

[c5, c6, c7, c8] = solve([f1 == 0, f2 == 0, f3 == 0, f4 == 1], [c5, c6, c7, c8])

t0 = 1; tf = 2; dt = 0.01;
t = t0:dt:tf;

x_dot = c5*t.^4 + c6*t.^3 + c7;
x = c5/5*t.^5 + c6/4*t.^4 + c7*t + c8;
figure(1)
plot(t, x_dot, 'LineWidth', 3)
xlabel('time (s)')
ylabel('Optimal x_d_o_t(t)')
figure(2)
plot(t, x, 'LineWidth', 3)
xlabel('time (s)')
ylabel('Optimal x(t)')