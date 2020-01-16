
function       [YR,ZR,UR,R]=gen_sample(rho,Y,Z,U,n)

% This function samples data for the population and creates the sample for each simulation
 
 R=rand(n,1)<=rho;  % sampling indicator

 YR=Y(R,:);         % sampled outcome
 ZR=Z(R,:);         % sampled features
 UR=U(R,:);         % sampled treatment
