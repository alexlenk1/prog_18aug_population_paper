function [lcb_boot,ucb_boot,se_boot]=se_boots_calc(Y,X,Z,Nboot)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% calculating the eicker-huber-white variance for coefficient on X (scalar)

W=[X,Z];
[N,KW]=size(W);
WW=W'*W;
iWW=inv(WW);

beta=iWW*(W'*Y);                % least squares coefficient
theta=beta(1,1);
eps=Y-W*beta;                   % residual
We=W.*(eps*ones(1,KW));
Var_ehw=iWW*(We'*We)*iWW;       % robust variance
se_ehw=sqrt(Var_ehw(1,1));      % standard error for first coefficient


Tboot=zeros(Nboot,1);
for iboot=1:Nboot,
    Boot=floor(N*rand(N,1)+1);
    UBoot=U(Boot,:);
    YBoot=Y(Boot,:);
    ZBoot=Z(Boot,:);
    WBoot=[XBoot,ZBoot];
    WWBoot=WBoot'*WBoot; 
    iWWBoot=inv(WWBoot);

    betaBoot=iWWBoot*(WBoot'*YBoot);                % least squares coefficient
    thetaBoot=betaBoot(1,1);
    epsBoot=YBoot-WBoot*betaBoot;                   % residual
    WeBoot=WBoot.*(epBoots*ones(1,KW));
    Var_ehwBoot=iWWBoot*(WeBoot'*WeBoot)*iWWBoot;       % robust variance
    se_ehwBoot=sqrt(Var_ehwBoot(1,1));      % standard error for first coefficient
    t=thetaBoot/se_ehwBoot;
    Tboot(iboot,1)=t;
end

Tboot=sort(Tboot);
t975=Tboot(round(Nboot*0.975),1);
t025=Tboot(round(Nboot*0.975),1);

ucb_boot=theta-t025*se_ehw;
lcb_boot=theta-t975*se_ehw;
se_boot=(ucb_boot-lcb_boot)/(2*1.96);


end