
function       [YR,XR,ZR,UR,R]=gen_sample(rho,Y,X,Z,U);

% generating sample data
 
R=rand(n,1)<=rho;  % sampling indicator

YR=Y(R,:);         % sampled outcome
XR=X(R,:);         % sampled adjusted treatment
ZR=Z(R,:);         % sampled features
UR=U(R,:);         % sampled treatment
