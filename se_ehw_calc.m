function se_ehw=se_ehw_calc(Y,X,Z)

% Calculates the Eicker-Huber-White variance for coefficient on X (scalar)

W=[X,Z];
[N,KW]=size(W);
WW=W'*W;
iWW=inv(WW);

beta=iWW*(W'*Y);                % least squares coefficient
eps=Y-W*beta;                   % residual
We=W.*(eps*ones(1,KW));
Var_ehw=iWW*(We'*We)*iWW;       % robust variance
se_ehw=sqrt(Var_ehw(1,1));      % standard error for first coefficient

end

