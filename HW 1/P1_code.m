clc
clear
close all

G1 = tf([143 2.6 26 0 0], [12 0.1 1])
G2 = tf([1 -0.1 -1], [-143 -2.6 -26 0 0])
p1 = pole(G1)
z1 = zero(G1)
p2 = pole(G2)
z2 = zero(G2)