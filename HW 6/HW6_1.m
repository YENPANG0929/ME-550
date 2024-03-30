clc
clear
close all

A = [0 1 ; -2 -3]; B = [0 ; 1];
Q = 0; R = 0.5; F = 0;
ti = 0; tf = pi/2;
x0 = [1 ; 0]; xf = [2 ; 0];

G_ = zeros(2);
options=odeset('RelTol',1e-10);
[t, G] = ode45(@(t,G) Gdiff(t,G,A,B,R),[0 tf-ti],G_,options);
figure(1)
plot(t, G,'LineWidth',3)
title('Gain G(t)')
xlabel('Time')
ylabel('Gain G')
set(gca, 'FontSize', 20)

Gf  = reshape(G(end,:),2,2);
uu = @(t) inv(R)*B'*expm(A'*(tf-t))*inv(Gf)*(xf - expm(A*(tf-ti))*x0);
u = zeros(size(t));
for i = 1:1:length(t)
    u(i) = uu(t(i));
end
figure(2)
plot(t, u,'LineWidth',3)
title('State u(t)')
xlabel('Time')
ylabel('Input u')
set(gca, 'FontSize', 20)

[t, x] = ode45(@(t,x) xdiff(t,x,A,B,uu),[0 tf-ti],x0,options);
figure(3)
plot(t, x,'LineWidth',3)
title('State x(t)')
xlabel('Time')
ylabel('State x')
set(gca, 'FontSize', 20)

function Gdiff = Gdiff(t,G,A,B,R)
G = reshape(G,2,2);

Gdiff = A*G + G*A' + B*inv(R)*B';
Gdiff = Gdiff(:);
end

function xdiff = xdiff(t,x,A,B,u)
u = u(t);

xdiff = A*x + B*u;
end


