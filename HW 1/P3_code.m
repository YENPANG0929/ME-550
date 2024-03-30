clc
clear
close all

A = [0 0 1 0; 0 0 0 1; -1/11 1/11 -0.1/11 0.1/11; 1/11 -1/11 0.1/11 -0.1/11];
B = [0; 0; 12/143; -1/143];
C = [0 1 0 0];
D = 0;
Q = C'*C;
R = 1;
[K,S,P] = lqr(A,B,Q,R);
[num, den] = ss2tf(A-B*K, B, C, D);
sys = tf(num, den);
damp(sys)
[wn, zeta, poles] = damp(sys);
Ts = (-log(0.02))/(zeta(1)*wn(1))
figure(1)
hold on
step(sys)
stepinfo(sys)
plot([0, 100], [1.02, 1.02], '--r', 'LineWidth', 1.5);
plot([0, 100], [0.98, 0.98], '--r', 'LineWidth', 1.5);
%% 
figure(2)
R = 1e-24;
[K,S,P] = lqr(A,B,Q,R);
[num, den] = ss2tf(A-B*K, B, C, D);
sys = tf(num, den);
step(sys)
stepinfo(sys)