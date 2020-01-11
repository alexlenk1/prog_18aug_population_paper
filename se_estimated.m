function se_est=se_estimated(Y,X,Z,rho);

% variance is Gamma^-1 Delta Gamma^-1

N=length(Y);                % number of units in sample
R=[X,ones(N,1),Z];          % matrix of regressors, X cause, 1 and Z attributes
beta=inv(R'*R)*(R'*Y);      % least squares estimates
eps=Y-R*beta;               % residuals

% calculation of Gamma (3 by 3 matrix)
G11=-mean(X.*X);
G12=-mean(X);
G13=-mean(X.*Z);
G22=-1;
G23=-mean(Z);
G33=-mean(Z.*Z);
Gamma=[G11,G12,G13;G12,G22,G23;G13,G23,G33];


% two estimates for Delta
Delta_V=zeros(3,3);   % estimate of proposed Delta estimator, Delta_strat in paper
Delta_ehw=zeros(3,3); % estimate of ehw Delta

% old version, based on matching

if true,
index=1:N;
index=index';

for i=1:N,
    Zi=Z(i,1);
    DZ=Z-Zi;
    DZ(i,1)=1000;
    if false,
    Mi=index(abs(DZ)<=min(abs(DZ))+0.0001,1);
    [S,L]=sort(rand(length(Mi),1));
    Mi=Mi(L,1);
    Mi=Mi(1,1);
    dxe=eps(i,1)*X(i,1)-eps(Mi,1)*X(Mi,1);
    de=eps(i,1)-eps(Mi,1);
    dze=eps(i,1)*Z(i,1)-eps(Mi,1)*Z(Mi,1);
    Delta_V=Delta_V+[dxe*dxe,dxe*de,dxe*dze;dxe*de,de*de,de*dze;dze*dxe,dze*de,dze*dze]/2;
    end,
    xe=eps(i,1)*X(i,1);
    ee=eps(i,1);
    ze=eps(i,1)*Z(i,1);
    Delta_ehw=Delta_ehw+[xe*xe,xe*ee,xe*ze;xe*ee,ee*ee,ee*ze;xe*ze,ze*ee,ze*ze];
    end
Delta_V=Delta_V/N;
Delta_ehw=Delta_ehw/N;

Var=inv(Gamma)*(rho*Delta_V+(1-rho)*Delta_ehw)*inv(Gamma);
se_est1=sqrt(Var(1,1)/N);
end       
%se_est1=0;

       
% new version, based on matching
Y1=Y(Z==1,:);
R1=R(Z==1,:);
Z1=Z(Z==1,1);
eps1=eps(Z==1,:);
X1=X(Z==1,:);
xe1=eps1.*X1;
ze1=eps1.*Z1;
e1=eps1;
e1=e1-mean(e1);
ze1=ze1-mean(ze1);
xe1=xe1-mean(xe1);


Y0=Y(Z==-1,:);
R0=R(Z==-1,:);
Z0=Z(Z==-1,1);
eps0=eps(Z==-1,:);
X0=X(Z==-1,:);
xe0=eps0.*X0;
ze0=eps0.*Z0;
e0=eps0;
e0=e0-mean(e0);
ze0=ze0-mean(ze0);
xe0=xe0-mean(xe0);


Delta_strat=zeros(3,3);
Delta_strat=Delta_strat+[xe1'*xe1,xe1'*e1,xe1'*ze1;xe1'*e1,e1'*e1,e1'*ze1;ze1'*xe1,e1'*ze1,ze1'*ze1];
Delta_strat=Delta_strat+[xe0'*xe0,xe0'*e0,xe0'*ze0;xe0'*e0,e0'*e0,e0'*ze0;ze0'*xe0,e0'*ze0,ze0'*ze0];
Delta_strat=Delta_strat/N;


Var=inv(Gamma)*(rho*Delta_strat+(1-rho)*Delta_ehw)*inv(Gamma);
se_est2=sqrt(Var(1,1)/N);
 
se_est=[se_est1;se_est2];
%pause