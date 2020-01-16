function se_boot=se_boots_calc(Y,X,Z,Nboot,N,se_ehw,theta)

% This function calculates improved bootstrapped standard errors.

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
        t=thetaBoot/se_ehwBoot;                         % bootstrapped t-statistic
        Tboot(iboot,1)=t;
    end
       
    % We calculate the improved critical values for the t-statistic. 
    Tboot=sort(Tboot);
    t975=Tboot(round(Nboot*0.975),1);                   %97.5 percentile of boostrap t-statistic distribution 
    t025=Tboot(round(Nboot*0.975),1);                   %2.5 percentile of boostrap t-statistic distribution 
    
    % We calculate the improved boundaries of the confidence interval under heteroskedasticity.
    ucb_boot=theta-t025*se_ehw;
    lcb_boot=theta-t975*se_ehw;
    
    % We find the corresponding standard error that is consistent with our imporved confidence interval. 
    % Since the standard error can uniquely be defined from the boundaries of the confidence interval only under the assumption of normality,
    % we divide by the critical values of the t-statistic under normality.
    se_boot=(ucb_boot-lcb_boot)/(2*1.96);
