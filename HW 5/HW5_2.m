clc
clear
close all

%% Problem 1

% Sys
A = [0 1; -2 0]; B = [0; 3]; C = 1; D = 0;
% Weight
Q = 0; R = 0.5; F = [0.5 0; 0 0];
% Initial
x0 = [0; 1];
ti = 0; tf = pi/2;

% Kalman gain
[tK, K] = Psolve(A, B, Q, R, F, ti, tf);
% Plot K(t)
figure(1)
hold on
plot(tK, K(:, 1), 'r', 'LineWidth', 3)
plot(tK, K(:, 2), 'r--', 'LineWidth', 3)
hold off
title('Kalman Gain K(t)')
xlabel('Time')
ylabel('Gain K')
legend('K_1', 'K_2')
set(gca, 'FontSize', 20)

% ODE45
options = odeset('RelTol', 1e-10);
[t, x] = ode45(@(t, x) xdiff(t, x, flag, A, B, tK, K),[0 (tf - ti)], x0, options);

% Plot x(t)
figure(2)
hold on
plot(t, x(:, 1), 'b', 'LineWidth', 3)
plot(t, x(:, 2), 'b--', 'LineWidth', 3)
hold off
title('State x(t)')
xlabel('Time')
ylabel('State x')
legend('x_1', 'x_2')
set(gca, 'FontSize', 20)

% Plot u(t)
[m, n] = size(x);
[mR, nR] = size(R);
Kt = interp1(tK, K, t);
u = zeros(m, mR);
for jj = 1:1:m
    u(jj) = -Kt(jj, :)*(x(jj, :))';
end
figure(3)
plot(t, u, 'r', 'LineWidth', 3)
title('Input u(t)')
xlabel('Time')
ylabel('Input u')
set(gca, 'FontSize', 20)





%%

% Similar but with Q matrix
% Sys
A_ = [0 1; -2 0]; B_ = [0; 3]; C_ = 1; D_ = 0;
% Weight
Q_ = 5; R_ = 0.5; F_ = [0.5 0; 0 0];
% Initial
x0_ = [0; 1];
ti_ = 0; tf_ = pi/2;

% Kalman gain
[tK_, K_] = Psolve(A_, B_, Q_, R_, F_, ti_, tf_);
% Plot K(t) vs K_Q(t)
figure(4)
hold on
plot(tK, K(:, 1), 'b', tK_, K_(:, 1), 'r', 'LineWidth', 3)
plot(tK, K(:, 2), 'b--', tK_, K_(:, 2), 'r--', 'LineWidth', 3)
hold off
title('Kalman Gain K(t) vs K_Q(t)')
xlabel('Time')
ylabel('Gain K')
legend('K_1', 'K_Q_1', 'K_2', 'K_Q_2')
set(gca, 'FontSize', 20)

% ODE45
options = odeset('RelTol', 1e-10);
[t_, x_] = ode45(@(t_, x_) xdiff(t_,x_,flag,A_,B_,tK_,K_),[0 (tf_ - ti_)], x0_, options);

% Plot x(t) vs x_Q(t)
figure(5)
hold on
plot(t, x(:, 1), 'b', t_, x_(:, 1), 'r', 'LineWidth', 3)
plot(t, x(:, 2), 'b--', t_, x_(:, 2), 'r--', 'LineWidth', 3)
hold off
title('State x(t) vs x_Q(t)')
xlabel('Time')
ylabel('State x')
legend('x_1', 'x_Q_1', 'x_2', 'x_Q_2')
set(gca, 'FontSize', 20)

% Plot u(t) vs u_Q(t)
[m_, n_] = size(x_);
[mR_, nR_] = size(R_);
Kt_ = interp1(tK_, K_, t_);
u_ = zeros(m_, mR_);
for jj_=1:1:m_
    u_(jj_) = -Kt_(jj_, :)*(x_(jj_, :))';
end
figure(6)
plot(t, u, 'b', t_, u_, 'r', 'LineWidth', 3)
title('Input u(t) vs u_Q(t)')
xlabel('Time')
ylabel('Input u')
legend('u', 'u_Q')
set(gca, 'FontSize', 20)

% Optimal control from HW 4
figure(7)
tt = linspace(0, pi/2, 100);
uu = -3*(0.0597*cos(sqrt(2)*tt) + 0.0455*sin(sqrt(2)*tt));
plot(tt, uu, 'b', t, u, 'r', 'LineWidth', 3)
title('Input u_H_W_4(t) vs u(t)')
xlabel('Time')
ylabel('Input u')
legend('u_H_W_4', 'u')
set(gca, 'FontSize', 20)



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Function to find the gain matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,K] = Psolve(A,B,Q,R,F,ti,tf)
%PSolve - Compute continuous-time solution to the Riccati Equation solved backward in time.
%
%   [t,K] = Psolve(A,B,Q,R,ti,tf) computes the gain matrix K
%   K = -R^(1) B' P  
%   by solving for the symmetric Riccati matrix P
%   from the Riccati equation
%       .                           
%       P = - PA  -A'P  -Q + PBR^(-1)B'P 
%   with final condition, i.e. P(tf) = F. 
%   usig Backward integration
%       .                           
%       P = + PA  +A'P  +Q - PEP 
%   where E = BR^(-1)B'
%   and  initial Condition P(0) = F
%   then needs to be reversed in time 
%    See also Pdiff.
%
[m,n] = size(A); NT = n*(n+1)/2; 
E = B*inv(R)*B';
options=odeset('RelTol',1e-10);

% Fiding final P in a vector form using F
Ptf = zeros(NT,1);
k = 1;
for i=1:n
    for j=i:n
        Ptf(k) = F(i,j);
        k = k+1;
    end
end

% Solving for P in a vector form PV

[t,PV]=ode45(@(t,p) PVdiff(t,p,flag,A,B,Q,R,E,F),[0 (tf-ti)],Ptf,options);
%
%
% PV is in vector form, each row corresponds to row in time t
% flip the PV vector 
PV = flipud(PV);
% redefine the time vector
t = flipud(t); 
t = -t +(tf)*ones(size(t));
%
%
% computing the gain matrix K(t) as row vector
%
[mP,nP] = size(PV);
K = zeros(mP,n);
%
for jj = 1:1:mP
    % find P matrix at time jj
    Pjj = zeros(n);
    for i=1:n    
        for j=i:n
            k = i*n - i*(i-1)/2 - (n-j);
            Pjj(i,j) = PV(jj,k);
            Pjj(j,i) = Pjj(i,j);
        end
    end
    % computing the feedback gain matrix
K(jj,:) = inv(R)*B'*Pjj;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Function to find the solution to Riccati Eq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PVdiff = PVdiff(t,p,flag,A,B,Q,R,E,F)
%
%    Called by PSolve to compute the Riccati matrix P.
%
%    [t,GP] = ode45('gramdiff',[ti tf],zeros(NT,1),[],A,B,Q,R,E).
%
%    See also PSOLVE.
%
%

% Finding the P matrix PM
[m,n] = size(A);
NT = n*(n+1)/2;
PM = zeros(n);
%
%
% Finding the P matrix PM
for i=1:n    
    for j=i:n
        k = i*n - i*(i-1)/2 - (n-j);
        PM(i,j) = p(k);
        PM(j,i) = PM(i,j);
    end
end
%
%
% computing P matrix derivative 
% Backward in time
PM_diff = A'*PM + PM*A -PM*E*PM +Q;
%
%
% computing P vector derivative
PVdiff = zeros(NT,1);
k = 1;
for i=1:n
    for j=i:n
        PVdiff(k) = PM_diff(i,j);
        k = k+1;
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Function to find the solution to system with optimal control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xdiff = xdiff(t,x,flag,A,B,tK,K)
%
% solving xdot = AX + Bu

Kt = interp1(tK,K,t);
ut = -Kt*x;
xdiff = A*x +B*ut;
end