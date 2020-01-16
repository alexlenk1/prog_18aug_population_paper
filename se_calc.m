function [se_ehw,se_desc,se_causal_sample,se_causal]=se_calc(Y,X,Z,beta,rho,N)

% This function outputs estimates for several variance estimators given in the paper (see p. 22-23) 

  H=X'*X; %denoted as capital Gamma in the paper 

  eps = Y-[X Z]*beta %calculating regression residuals (as on p.22) 
  Xeps=X.*eps;
  G=(inv(Z'*Z)*(Z'*Xeps));

  Delta_ehw=(Xeps'*Xeps)/N;
  Delta_z=(Xeps-Z*G)'*(Xeps-Z*G)/N;

  % These are estimated sqrt(N) variances
  V_ehw=inv(H)*Delta_ehw*inv(H);
  V_desc = (1-rho)*V_ehw;
  V_causal_sample=inv(H)*Delta_z*inv(H);
  V_causal=rho*V_causal_sample + (1-rho)*V_ehw;


  % We now calculate the estimated variances of the estimators (no longer sqrt(N))
  se_ehw=sqrt(V_ehw/N);
  se_desc=sqrt(V_desc/N);
  se_causal=sqrt(V_causal/N);
  se_causal_sample = sqrt(V_causal_sample/N);
