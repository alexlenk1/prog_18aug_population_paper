% generating population data

function [Y,U,Z,theta_causal,theta_desc,W,Omega]=gen_pop(M,tau,psi,beta,gamma,xi);



Z=[ones(n,1),randn(n,1);
ZZ=Z'*Z/n;
Mean_Z=mean(Z);

eta=randn(n,1);
U=Z*xi+randn(n,1);
X=U-Z*Xi;
[nZ,KZ]=size(Z);
Y=U*theta+Z*gamma+((U*ones(1,KZ)).*Z)*psi+eta;

W=[Y X Z]'*[Y X Z]/n;
EXX=1;     % expected value of X*X
EXZ=[0 0]; % expected value of X*Z
EXY=theta+Mean_Z*phi;
EYZ=xi'*ZZ*theta+ZZ*gamma+Mean_Z*xi'*ZZ*psi+Z'*eta/n;
EYY=;
Omega=[EYY EXY EYZ; EXY EXX EXZ; EYZ' EXZ' Z'*Z/n];

Exp_U=mean_Z*xi;

% the estimands
theta_causal=tau+(exp_U*ones(1,KZ))*(Mean_Z')*gamma;
theta_desc=inv([U Z]'*[U Z])*([U Z]'*Y);
theta_desc=theta_desc(1,1);
