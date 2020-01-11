% generating sample data
function       [YR,XR,ZR,W,UR,R,TOmega,theta_desc,theta_causal_sample]=gen_sample(alpha,phi,Z,rho,xi);
  

n=length(alpha);   % population size

U=Z*xi+randn(n,1); % -1+2*(rand(n,1)<0.5); % raw treatment
X=U-Z*xi;          % adjusted treatment


Y=alpha+X.*phi;    % realized outcome
R=rand(n,1)<=rho;  % sampling indicator
N=sum(R);          % sample size

YR=Y(R,:);         % sampled outcome
XR=X(R,:);         % sampled adjusted treatment
ZR=Z(R,:);         % sampled features
UR=U(R,:);         % sampled treatment

W=[YR XR ZR]'*[YR XR ZR]/N;

beta=inv([X Z]'*[X Z])*([X Z]'*Y);
theta_desc=beta(1,1);

theta_causal_sample=mean(phi(R,1));

TOmega=0; %[EYYR EYXR EYZR; EYXR EXXR EXZR; EYZR' EXZR' ZZR];