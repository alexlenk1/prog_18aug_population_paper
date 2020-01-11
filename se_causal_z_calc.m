function se_causal_z=se_causal_z_calc(Y,U,Z,rho)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% calculating the eicker-huber-white variance for coefficient on X (scalar)

beta=inv([U Z]'*[U Z])*([U Z]'*Y);
eps=Y-[U Z]*beta;
Lambda=(inv(Z'*Z)*(Z'*U));
X=U-Z*Lambda;
H=X'*X;
Xeps=X.*eps;

G=(inv(Z'*Z)*(Z'*Xeps));
rsq=(Z*G)'*(Z*G)/(Xeps'*Xeps);

Delta_z=(Xeps-Z*G)'*(Xeps-Z*G);
Delta_ehw=Xeps'*Xeps;
Delta_causal_z=rho*Delta_z+(1-rho)*Delta_ehw;
V_ehw=inv(H)*Delta_ehw*inv(H);
se_ehw=sqrt(V_ehw);
V_z=inv(H)*Delta_z*inv(H);
V_causal_z=inv(H)*Delta_causal_z*inv(H);


se_causal_z=sqrt(V_causal_z);

end

