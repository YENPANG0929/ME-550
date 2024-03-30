clc
clear
close all

A = [0 0 1 0; 0 0 0 1; -1/11 1/11 -0.1/11 0.1/11; 1/11 -1/11 0.1/11 -0.1/11];
B = [0; 0; 12/143; -1/143];
C = [0 1 0 0];
D = 0;
[num, den] = ss2tf(A,B,C,D);
sys = tf(num, den)
poles = pole(sys)
zeros = zero(sys)