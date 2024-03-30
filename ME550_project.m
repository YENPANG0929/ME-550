clc
clear
close all



% Time
num_agents = 12;
t0 = 0; tf= 10; dt = 0.01;
t = t0:dt:tf;

% Init
x_range = [-4, 4]; % x range of the environment
y_range = [-4, 4]; % y range of the environment

initial_x = (x_range(2) - x_range(1)) * rand(1, num_agents) + x_range(1); 
initial_y = (y_range(2) - y_range(1)) * rand(1, num_agents) + y_range(1);

x0_all = [initial_x; initial_y];
xf = [0; 0];

% System
A = [0 2; -3 -1];
B = [0; 2];
C = [1 0; 0 0];
Q = 0.5*eye(2); % change this to make tracking better
r = 1;
R = 0.5*r;
F = 10*Q;

% 
for a = 1:num_agents
    x0 = x0_all(:, a);
    % Simplified pri with Gaussian noises instead
    noise_std_dev = 0.1; % Standard deviation for noises
    pri = @(t) [0.1; 0.07]*t + [cos(0.02*t) -sin(0.02*t); cos(0.02*t) cos(0.02*t)] * normrnd(0, noise_std_dev, 2, 1);

    [tK,K,G] = Psolve(A,B,C,Q,R,F,t0,tf,pri); 
    options=odeset('RelTol',1e-10);
    [tx,x]=ode45(@(tx,x) xdiff(tx,x,flag,A,B,R,tK,K,G),[t0 tf],x0,options);
    [m,n]=size(x); % Length of time, x is m
    [mR,nR] = size(R); % number of inputs = mR
    Kt = interp1(tK,K,tx); % interpolate the K matrix
    Gt = interp1(tK,G,tx); % interpolate the G matrix
    
    u = zeros(m,mR); % initialize the input
    for jj=1:1:m
    u(jj) = -Kt(jj,:)*(x(jj,:))' + inv(R)*B'*(Gt(jj,:))'; % 
    end
    
    % Formation vectors
    % Circle
    hx = cos((2*(a - 1)*pi/num_agents) + 2*tx);
    hy = sin((2*(a - 1)*pi/num_agents) + 2*tx);
    
    % Star
%     hx = cos(a * tx) .* cos(tx);
%     hy = cos(a * tx) .* sin(tx);

    % Lemniscate of Gerono
%     hx = a * cos(tx);
%     hy = a * sin(2 * tx);

    % Heart
%     hx = a * (16 * sin(tx).^3);
%     hy = a * (13 * cos(tx) - 5 * cos(2 * tx) - 2 * cos(3 * tx) - cos(4 * tx));

    % Local formation references
    xri = x(:,1) - hx; 
    yri = x(:,2) - hy;

    figure(1)
    grid on
    plot3(x(:,1), x(:,2), tx, 'LineWidth', 3)
    hold on
    scatter3(x(1, 1), x(1, 2), tx(1), 100, 'MarkerEdgeColor', 'k', 'MarkerFaceColor','b')
    scatter3(x(end, 1), x(end, 2), tx(end), 100, 'MarkerEdgeColor', 'k', 'MarkerFaceColor','r')
    set(gca,'FontSize',20)

    figure(2)
    grid on
    plot3(xri, yri, tx, 'LineWidth', 3)
    hold on
    scatter3(xri(1), yri(1), tx(1), 100, 'MarkerEdgeColor', 'k', 'MarkerFaceColor','b')
    scatter3(xri(end), yri(end),tx(end), 100, 'MarkerEdgeColor', 'k', 'MarkerFaceColor','r')
    set(gca,'FontSize',20)

    figure(3)
    grid on
    hold on
    scatter(xri(1), yri(1), 100, 'filled', 'MarkerEdgeColor', 'k')
    set(gca,'FontSize',20)

    figure(4)
    grid on
    hold on
    scatter(xri(round(2*end/5)), yri(round(2*end/5)), 100, 'filled', 'MarkerEdgeColor', 'k')
    set(gca,'FontSize',20)
    
    figure(5)
    grid on
    hold on
    scatter(xri(round(3*end/5)), yri(round(3*end/5)), 100, 'filled', 'MarkerEdgeColor', 'k')
    set(gca,'FontSize',20)

    figure(6)
    grid on
    hold on
    scatter(xri(round(4*end/5)), yri(round(4*end/5)), 100, 'filled', 'MarkerEdgeColor', 'k')
    set(gca,'FontSize',20)

    figure(7)
    grid on
    hold on
    scatter(xri(round(5*end/5)), yri(round(5*end/5)), 100, 'filled', 'MarkerEdgeColor', 'k')
    set(gca,'FontSize',20)

end

figure(1)
title('The position trajectories without formation vector')
xlabel('q_i_,_X(t)')
ylabel('q_i_,_Y(t)')
zlabel('t(s)')
legend('Agent1', '', '', 'Agent2', '', '', 'Agent3', '', '', 'Agent4', '', '', 'Agent5', '', '', 'Agent6', '', '', 'Agent7', '', '', 'Agent8', '', '', 'Agent9', '', '', 'Agent10', '', '', 'Agent11', '', '', 'Agent12', 'Location', 'Best')

figure(2)
title('The position trajectories with formation vector')
xlabel('q_i_,_X(t)')
ylabel('q_i_,_Y(t)')
zlabel('t(s)')
legend('Agent1', '', '', 'Agent2', '', '', 'Agent3', '', '', 'Agent4', '', '', 'Agent5', '', '', 'Agent6', '', '', 'Agent7', '', '', 'Agent8', '', '', 'Agent9', '', '', 'Agent10', '', '', 'Agent11', '', '', 'Agent12', 'Location', 'Best')

figure(3)
title('Snapshot at t = 0')
xlabel('q_i_,_X(t)')
ylabel('q_i_,_Y(t)')
zlabel('t(s)')
legend('Agent1', 'Agent2', 'Agent3', 'Agent4', 'Agent5', 'Agent6', 'Agent7', 'Agent8', 'Agent9', 'Agent10', 'Agent11', 'Agent12', 'Location', 'Best')

figure(4)
title('Snapshot at t = 2.5')
xlabel('q_i_,_X(t)')
ylabel('q_i_,_Y(t)')
zlabel('t(s)')
legend('Agent1', 'Agent2', 'Agent3', 'Agent4', 'Agent5', 'Agent6', 'Agent7', 'Agent8', 'Agent9', 'Agent10', 'Agent11', 'Agent12', 'Location', 'Best')

figure(5)
title('Snapshot at t = 5')
xlabel('q_i_,_X(t)')
ylabel('q_i_,_Y(t)')
zlabel('t(s)')
legend('Agent1', 'Agent2', 'Agent3', 'Agent4', 'Agent5', 'Agent6', 'Agent7', 'Agent8', 'Agent9', 'Agent10', 'Agent11', 'Agent12', 'Location', 'Best')

figure(6)
title('Snapshot at t = 7.5')
xlabel('q_i_,_X(t)')
ylabel('q_i_,_Y(t)')
zlabel('t(s)')
legend('Agent1', 'Agent2', 'Agent3', 'Agent4', 'Agent5', 'Agent6', 'Agent7', 'Agent8', 'Agent9', 'Agent10', 'Agent11', 'Agent12', 'Location', 'Best')

figure(7)
title('Snapshot at t = 10')
xlabel('q_i_,_X(t)')
ylabel('q_i_,_Y(t)')
zlabel('t(s)')
legend('Agent1', 'Agent2', 'Agent3', 'Agent4', 'Agent5', 'Agent6', 'Agent7', 'Agent8', 'Agent9', 'Agent10', 'Agent11', 'Agent12', 'Location', 'Best')



%% Helper Functions
function [tK,K,G] = Psolve(A,B,C,Q,R,F,ti,tf,z)
[n,m] = size(B);
NT = n*(n+1)/2;
E = B*inv(R)*B';
options=odeset('RelTol',1e-12);
% Fiding final P in a vector form using F
Ptf = zeros(NT,1);
PMtf = C'*F*C;
k = 1;
for i=1:n
for j=i:n
Ptf(k) = PMtf(i,j);
k = k+1;
end
end
gtf = [0;0];
% combined vector
PGtf = [Ptf; gtf];
% Solving for P in a vector form PV
[t,PVG]=ode45(@(t,p) PVGdiff(t,p,flag,A,B,C,Q,R,F,tf,z),[tf ti],PGtf,options);
%
%
% PV is in vector form, each row corresponds to row in time t
% flip the PV vector
PVG = flipud(PVG);
% redefine the time vector
t = flipud(t);
% t = -t +(tf)*ones(size(t));
%
%
% computing the gain matrix K(t) as row vector
%
PV = PVG(:,1:NT);
G = PVG(:,NT+1:NT+n);
% plot the riccati matrix terms
% nfig = nfig+1; figure(nfig)
figure(10)
grid on
hold on
subplot(1,2,1);
plot(t,PV,'LineWidth',3)
title('Riccati matrix P(t)')
xlabel('time')
ylabel('P(t)')
legend('P_{11}','P_{12}','P_{22}','Location','best')
set(gca,'FontSize',20)
% plot the riccati matrix terms
% nfig = nfig+1; figure(nfig)
subplot(1,2,2);
plot(t,G,'LineWidth',3)
title('Riccati matrix g(t)')
xlabel('time')
ylabel('g(t)')
legend('g_{1}','g_{2}','Location','best')
set(gca,'FontSize',20)
[mP,nP] = size(PVG);
K = zeros(mP,n);
G = zeros(mP,n);
tK =t;
%
for jj = 1:1:mP
% find P matrix at time jj
Pjj = zeros(n);
for i=1:n
for j=i:n
k = i*n - i*(i-1)/2 - (n-j);
Pjj(i,j) = PVG(jj,k);
Pjj(j,i) = Pjj(i,j);
end
end
% computing the feedback gain matrix
K(jj,:) = inv(R)*B'*Pjj;
G(jj,:) = PVG(jj,NT+1:NT+n);
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Function to find the solution to Riccati Eq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PVGdiff = PVGdiff(t,p,flag,A,B,C,Q,R,F,tf,z)
%
% Called by PSolve to compute the Riccati matrix P.
%
% See also PSOLVE.
%
%
% Finding the P matrix PM
[m,n] = size(A);
NT = n*(n+1)/2;
PM = zeros(n);
E = B*inv(R)*B';
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
PM_diff = -(A'*PM + PM*A -PM*E*PM +C'*Q*C);
%
zt = z(t);
gt = p(NT+1: NT+n,1);
%
Gdiff = -((A -E*PM)'*gt +C'*Q*zt);
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
PVGdiff = [PVdiff; Gdiff];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Function to find the solution to system with optimal control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xdiff = xdiff(t,x,flag,A,B,R,tK,K,G)
%
% solving xdot = AX + Bu
Kt = interp1(tK,K,t);
Gt = interp1(tK,G,t);
ut = -Kt*x +inv(R)*(B')*Gt';
xdiff = A*x +B*ut;
end
