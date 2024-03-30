clc
clear
close all

nfig = 0;
%% Problem 1
% describe the system and PI matrices
    % System matrices 
    A = @(k) [0.1 (1+0.1*k) 0; 0 0.2 0.1; 0 0 0.2]; B=[1;2;10]; C = [1 0 0]; D =0;
    % cost function matrices
    Q = 100*[1]; R = 1; F = 100*[1];
    % initial and final time  
    k0 = 1; kf =10;
    % initial conditions 
    x0 = [1;2;10]; 
    

tkf=k0:1:(kf); 
%z=2*ones(size(tkf));
z= @(k) 2*(k); % another trajectory
%z=sin(tkf); % another trajectory

%return
% Find the controller gain K and g 
[Lx, Lg, G] = Psolve(A,B,C,Q,R,F,k0,kf,tkf,z);


tk=k0:1:(kf-1); 
% plot the control gain K
nfig = nfig+1; figure(nfig)
plot(tk,Lx,'ko',tk,Lx,'k:','LineWidth',2)
xlabel('time step k')
ylabel('Kalman gain L_x')
set(gca,'FontSize',20)

% plot the control gain K
nfig = nfig+1; figure(nfig)
plot(tk,Lg,'ko',tk,Lg,'k:','LineWidth',2)
xlabel('time step k')
ylabel('L_g')
set(gca,'FontSize',20)

% plot the control gain K
nfig = nfig+1; figure(nfig)
plot(tkf,G,'ko',tkf,G,'k:','LineWidth',2)
xlabel('time step k')
ylabel('G')
set(gca,'FontSize',20)

%return

% simulate the sytem response 
[X,U] = Syssolve(A,B,Lx,Lg,G,x0,k0,kf);
tx=k0:1:kf; 

% plot the state x
nfig = nfig+1; figure(nfig)
plot(tx,X(1,:),'bo',tx,X(1,:),'b:',tx,X(2,:),'go',tx,X(2,:),'g:','LineWidth',2)
xlabel('time step k')
ylabel('State x')
set(gca,'FontSize',20)


nfig = nfig+1; figure(nfig); clf
plot(tk,U,'ro',tk,U,'r:','LineWidth',2)
xlabel('time step k')
ylabel('Input u')
set(gca,'FontSize',20)















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Function to find the gain matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Lx,Lg,G] = Psolve(A,B,C,Q,R,F,k0,kf,tkf,z)
%
%
% initialize the Riccati matrices
[m,n] = size(A(kf));
[mz,nz] = size(z);
Pk1 = C'*F*C;
gk1 = C'*F*z(kf);
Pk   = zeros(n);

E = B*inv(R)*B';
V = C'*Q*C; 
W = C'*Q; 

Lx = [];
Lg = []; 
G = [gk1'];
II = eye(n); 
JL = 1:1:(kf-k0);
    for jj=length(JL):-1:1;
        Ak = A(JL(jj));
        zk = z(JL(jj));
        Pk = Ak'*Pk1*(inv(II+E*Pk1))*Ak + V; 
        lx = inv(R +B'*Pk1*B)*B'*Pk1*Ak;
        lg = inv(R +B'*Pk1*B)*B';
        g = (Ak'-Ak'*Pk1*(inv(II+E*Pk1))*E)*gk1 +W*zk;
        Pk1 = Pk;
        gk1 = g;
        Lx = [lx; Lx]; % stack the gain matrices
        Lg = [lg; Lg]; % stack the gain matrices
        G = [g'; G]; % stack the gain matrices
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Function to find the gain matrix
%%%% [X,U] = Syssolve(A,B,Lx,Lg,g,x0,k0,kf);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,U] = Syssolve(A,B,Lx,Lg,G,x0,k0,kf)
%
%
[m,n] = size(A(kf));
x = x0;
X = [x];
U=[];
    for jj=1:1:(kf-k0)
        u = -Lx(jj,:)*x +Lg(jj,:)*G(jj+1,:)';
        xp1 = A(jj) *x +B*u;
        X = [X xp1]; % stack the state vectors
        U = [U u];
    end
end