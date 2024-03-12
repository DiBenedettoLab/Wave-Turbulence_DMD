function [Phi,omega,lambda,b,Xdmd] = DMDselective(X1,X2,r,dt)
% function [Phi,omega,lambda,b,Xdmd] = DMD(X1,X2,r,dt)
% Computes the DynamicModeDecomposition of X1, X2
% INPUTS:
% X1 = X, data matrix
% X2 = X , shifted data matrix
% Columns of X1 and X2 are state snapshots
% r = target rank of SVD r(1) start, r(2) finish
% dt = time stepadvancing X1 to X2 (X to X’)
% OUTPUTS:
% Phi, the DMD modes
% omega, thecontinuous-time DMDeigenvalues 
% lambda, the discrete -time DMDeigenvalues
% b, a vector of magnitudes of modes Phi
% Xdmd, the datamatrix reconstructed by Phi,omega, b
%% DMD
[U,S,V]=svd(X1,'econ');
r(2)=min(r(2),size(U,2));U_r=U(:, r(1):r(2)); %truncate to rank-r

S_r=S(r(1):r(2),r(1):r(2));
V_r=V(:, r(1):r(2));
% Project A into dominant singular vectors
Atilde=U_r' *X2*V_r/S_r;% low-rankdynamics I added the transpose to U_r
[W_r,D]=eig(Atilde);
Phi=X2*(V_r/S_r)*W_r;% DMD modes maybe delete parenthesis

lambda=diag(D);% discrete -time eigenvalues
omega=log(lambda)/dt;% continuous-time eigenvalues (for predictions not in dt)

%% Compute DMDmodeamplitudes b
x1=X1(:, 1); % initial values
b=Phi\x1;

% another way of computing b
alpha1 = S_r * V_r(1,:)';
b = (W_r * D)\alpha1;

%% DMD reconstruction
mm1=size(X1,2);% mm1 = m - 1
time_dynamics=zeros(r(2)-r(1)+1,mm1);
t= (0:mm1-1)*dt;% time vector 
for iter = 1:mm1
    time_dynamics(:,iter) = (b.*exp(omega*t(iter)));
end
Xdmd=Phi*time_dynamics;