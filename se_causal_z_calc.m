function se_causal_z=se_causal_z_calc(Y,X,Z,beta_tilde,rho,N)

% This function outputs estimates for several variance estimators given in the paper 

% 1.Calculating the Eicker-Huber-White variance for coefficient on X (scalar)

%eps=Y-[U Z]*beta; NOT CONSISTENT WITH definition on p.22

H=X'*X; %denoted as capital Gamma in the paper 
eps = Y-[X Z]*beta_tilde %calculating regression residulas 
Xeps=X.*eps;
G=(inv(Z'*Z)*(Z'*Xeps));

Delta_ehw=(Xeps'*Xeps)/N;
Delta_z=(Xeps-Z*G)'*(Xeps-Z*G)/N;
Delta_causal_z=rho*Delta_z+(1-rho)*Delta_ehw;


V_ehw=inv(H)*Delta_ehw*inv(H); % This is the estimated sqrt(N) variance
V_desc = (1-rho)*V_ehw;
V_causal_sample=inv(H)*Delta_z*inv(H);
V_causal=inv(H)*Delta_causal_z*inv(H);

se_ehw=sqrt(V_ehw);
se_hat_desc=sqrt(V_desc)
se_causal_z=sqrt(V_causal_z);
se_causal_sample_z = sqrt(



V_z=inv(H)*Delta_z*inv(H);
V_causal_z=inv(H)*Delta_causal_z*inv(H);


se_causal_z=sqrt(V_causal_z);

end

