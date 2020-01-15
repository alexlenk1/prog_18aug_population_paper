function [lcb_boot,ucb_boot,se_boot]=se_boots_calc(Y,X,Z,Nboot,N,se_ehw,theta)

% This function calculates bootstrapped standard errors.

    Tboot=zeros(Nboot,1);
    for iboot=1:Nboot,
        Boot=floor(N*rand(N,1)+1);
        XBoot=X(Boot,:);
        YBoot=Y(Boot,:);
        ZBoot=Z(Boot,:);
        WBoot=[XBoot,ZBoot];
        WWBoot=WBoot'*WBoot; 
        iWWBoot=inv(WWBoot);

        betaBoot=iWWBoot*(WBoot'*YBoot);                % least squares coefficient
        thetaBoot=betaBoot(1,1);
        epsBoot=YBoot-WBoot*betaBoot;                   % residual
        WeBoot=WBoot.*(epBoots*ones(1,KW));
        Var_ehwBoot=iWWBoot*(WeBoot'*WeBoot)*iWWBoot;   % robust variance
        se_ehwBoot=sqrt(Var_ehwBoot(1,1));              % standard error for thetaboot coefficient
        t=thetaBoot/se_ehwBoot;                         % t-statistic
        Tboot(iboot,1)=t;
    end

    Tboot=sort(Tboot);
    t975=Tboot(round(Nboot*0.975),1);                   %97.5 percentile of boostrap t-statistic distribution 
    t025=Tboot(round(Nboot*0.975),1);                   %2.5 percentile of boostrap t-statistic distribution 

    ucb_boot=theta-t025*se_ehw;
    lcb_boot=theta-t975*se_ehw;
    se_boot=(ucb_boot-lcb_boot)/(2*1.96);


    end
