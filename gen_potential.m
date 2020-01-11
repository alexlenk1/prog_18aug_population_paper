% generating population data

function [Z,alpha,phi]=gen_potential(n,psi,gamma,xi,sig_eps,sig_eta,KZ);

% generate the potential outcome functions
% and covariates Z
% Y(u)=alpha+u*phi
% alpha is n-vector
% phi is n-vector
% n is population size


Z=[ones(n,1),randn(n,KZ)];   % exogenous fixed features
ZZ=Z'*Z/n; 
Mean_Z=mean(Z);

eta=randn(n,1)*sig_eta;
eps=randn(n,1)*sig_eps;

alpha=Z*gamma+eta;  % intercept for potential outcomes, depends on features
phi=Z*psi+eps;      % treatment effect, depents on features



if false,
U=Z*xi+randn(n,1);       % treatment
X=U-Z*xi;                % adjusted treatment
[nZ,KZ]=size(Z);
Y=alpha+U*phi;  % realized outcome


W=[Y X Z]'*[Y X Z]/n;
EXX=1;     % expected value of X*X
EXZ=[0 0]; % expected value of X*Z
EXY=theta+Mean_Z*phi;
EYZ=xi'*ZZ*theta+ZZ*gamma+Mean_Z*xi'*ZZ*psi+Z'*eta/n;
EYY=0;
Omega=[EYY EXY EYZ; EXY EXX EXZ; EYZ' EXZ' Z'*Z/n];

Exp_U=mean_Z*xi;


% the estimands
theta_causal=tau+(exp_U*ones(1,KZ))*(Mean_Z')*gamma;
theta_desc=inv([U Z]'*[U Z])*([U Z]'*Y);
theta_desc=theta_desc(1,1);
end