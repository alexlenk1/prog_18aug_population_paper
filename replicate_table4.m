% This program reproduces the Simulation Table 4 on p.39
% The functions that are called are:

% 1) gen_population.m
% Generates population of size n

% 2) gen_sample.m
% Generates sample of expected size (rho*N) for each simulation

% 3) se_calc.m
% Calculates estimates for different standard error estimators 

% 4) se_boots_calc.m
% Calculates improved bootstrapped standard erorrs 

% 5) gtable.m
% Creates table


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% IMPLEMENTATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear   % ensure that the workspace is empty

load seed_oct2 % loading random seed for reproducibility

Nboot=1000;           % number of bootstrap replications

% parameters for simulations
Nsim=50000;      % number of simulations
csim=1000;       % number of simulations before output is sent to screen for monitoring
EN=1000;         % expected number of observations in the sample.
                 % if the population size is n, and the sampling rate is
                 % rho, EN is approximately rho*n. EN and rho are fixed in the
                 % simulations, and the population size is calculated as
                 % EN/rho.
KZ=2;            % number of fixed characteristics Z, beyond intercept.

tabel=zeros(39,3);   % matrix containing output for tables
                     % there are seven simulation designs

for simulation_design=1:7
    if simulation_design==1
       KZ=1;      
       EN=1000;   
       rho=0.01; 
       n=round(EN/rho);
       psi=zeros(KZ+1,1);   
       psi(2,1)=2;
       gamma=zeros(KZ+1,1);
       lambda=zeros(KZ+1,1);
       sig_xi=1;
       sig_theta=1;  
    end
    
    if simulation_design==2
       KZ=10;      
       EN=1000;   
       rho=0.01; 
       n=round(EN/rho);
       psi=zeros(KZ+1,1);   
       psi(2,1)=2;
       gamma=zeros(KZ+1,1);
       lambda=zeros(KZ+1,1);
       sig_xi=1;
       sig_theta=1;  
    end

    if simulation_design==3
       KZ=1;     
       EN=100;  
       rho=0.01;  
       n=round(EN/rho);
       psi=zeros(KZ+1,1);   
       psi(2,1)=2;
       gamma=zeros(KZ+1,1);
       lambda=zeros(KZ+1,1);
       sig_xi=1;
       sig_theta=1;  
    end
    
    if simulation_design==4
       KZ=1;     
       EN=1000;  
       rho=1; 
       n=round(EN/rho);
       psi=zeros(KZ+1,1);   
       psi(2,1)=2;
       gamma=zeros(KZ+1,1);
       lambda=zeros(KZ+1,1);
       sig_xi=1;
       sig_theta=1; 
    end
    
    if simulation_design==5
       KZ=1;      
       EN=1000;   
       rho=0.01; 
       n=round(EN/rho);
       psi=zeros(KZ+1,1);   
       gamma=zeros(KZ+1,1);
       lambda=zeros(KZ+1,1);
       sig_xi=1;
       sig_theta=1; 
    end
    
    if simulation_design==6
       KZ=1;     
       EN=1000;  
       rho=0.01;  
       n=round(EN/rho);
       psi=zeros(KZ+1,1);   
       psi(2,1)=2;
       gamma=zeros(KZ+1,1);
       lambda=zeros(KZ+1,1);
       sig_xi=1;
       sig_theta=0;  
    end
    
    if simulation_design==7
       KZ=1;    
       EN=1000;   
       rho=0.01; 
       n=round(EN/rho);
       psi=zeros(KZ+1,1);   
       gamma=zeros(KZ+1,1);
       lambda=zeros(KZ+1,1);
       sig_xi=1;
       sig_theta=0;  
    end
    
    % Given the simulation design, we first generate the population of
    % units. There are n units in ths population, each characterized by
    % three variables:
    % 1. Z_i, a vector of length KZ+1, with the first element equal to 1
    % 2. xi_i, the intercept in the potential outcome function 
    % 3. theta_i, the slope coefficient in the potential outcome function
   
    [Z,xi,theta]=gen_population(n,psi,gamma,lambda,sig_theta,sig_xi,KZ);
    
    % Given our population, we can calculate the corresponding realized outcomes for each unit.
    % The potential outcome function is:
    % y(u)=xi_i+theta_i*u
    % After removing the correlation between U and X, the potential outcome function is: 
    % y(x)=xi_i+theta_i*x
    
    U=Z*lambda+randn(n,1); % raw treatment, distributed N(Z*lambda, 1)
    X=U-Z*lambda;          % adjusted treatment (after removing correlation from Z); in our case since lambda=0, X=U

    Y=xi+X.*theta;         % realized outcome
        
    % We now calculate the estimands for theta causal and theta descriptive: 
    
    theta_causal=mean(theta); % [This is implied by formula (3.5) since X,Z standard normal and uncorrelated and gamma=0]

    beta=inv([X Z]'*[X Z])*([X Z]'*Y); % Formula (3.3)
    theta_desc=beta(1,1);    
      
    % We now calculate the true variances of the different estimators.
    % Notice these are sqrt(N) variances as on p.18.
    Delta_ehw=1+3*(psi'*psi+sig_theta*sig_theta);
    Delta_Z=Delta_ehw-(psi'*psi);
    Delta_cond=Delta_ehw-(psi'*psi+sig_theta*sig_theta);
   
    H=1;  % variance of X, so with X being standard normal, this is equal to 1 (in paper, H is denoted as capital Gamma)
    V_ehw=inv(H)*Delta_ehw*inv(H);
    V_desc=(1-rho)*inv(H)*Delta_ehw*inv(H);
    V_causal_sample=inv(H)*Delta_cond*inv(H);
    V_causal=inv(H)*(rho*Delta_cond+(1-rho)*Delta_ehw)*inv(H);
    V_Z_causal_sample=inv(H)*Delta_Z*inv(H);
    V_Z_causal=inv(H)*(rho*Delta_Z+(1-rho)*Delta_ehw)*inv(H);
    
    % We now calculate the (expected) true standard errors of the estimators (these are no longer sqrt(N)).
    % These are expected true standard errors as we divide by expected rather than true sample size. 
    se_ehw=sqrt(V_ehw/(EN));
    se_desc=sqrt(V_desc/(EN));
    se_causal_sample=sqrt(V_causal_sample/(EN));
    se_causal=sqrt(V_causal/(EN));
    se_Z_causal_sample=sqrt(V_Z_causal_sample/(EN));
    se_Z_causal=sqrt(V_Z_causal/(EN));
    
    %Define Matrices for Storing Output from Each Simulation 
    
    theta_out=zeros(Nsim,4);  % this records estimated theta and how it deviates from the 3 theta estimators defined above
    out=zeros(Nsim,11); %  this records different variance estimators (both the true ones and the estimated ones)
    coverage=zeros(Nsim,33); % this records coverage rates for all variance-estimator combinations
    
    for isim=1:Nsim
        [YR,ZR,UR,R]=gen_sample(rho,Y,Z,U,n);
        
        N = sum(R); %sample size
        rho_hat = N/n;
        
        %Calculating true theta causal sample 
        theta_causal_sample=mean(theta(R,1)); % This is implied by formula (3.5) since X,Z standard normal and uncorrelated and gamma=0
        
        % Calculating estimated X (which is U net of correlation with X, see p.11)
        Lambda_hat=(inv(ZR'*ZR)*(ZR'*UR));
        XR=UR-ZR*Lambda_hat;
        
        % Calculating theta hat and gamma hat as on p.11
        beta_hat=inv([XR ZR]'*[XR ZR])*([XR ZR]'*YR); 
        hat_theta=beta(1,1);    % Remember that theta_hat = theta_tilde (where theta_tilde obtained in a regression of YR on UR and ZR)
        
        % Sample standard errors
        [se_hat_ehw,se_hat_desc,se_hat_causal_sample,se_hat_causal]=se_calc(YR,XR,ZR,beta_hat,rho_hat,N);
        se_hat_boot=se_boots_calc(YR,XR,ZR,Nboot,se_hat_ehw,hat_theta);
        
        % Calculating nominal coverage (ie, using (expected) true variances) 
        in_ehw_desc=abs(hat_theta-theta_desc)<1.96*se_ehw;
        in_ehw_causal_sample=abs(hat_theta-theta_causal_sample)<1.96*se_ehw;
        in_ehw_causal=abs(hat_theta-theta_causal)<1.96*se_ehw;
        
        in_desc_desc=abs(hat_theta-theta_desc)<1.96*se_desc;
        in_desc_causal_sample=abs(hat_theta-theta_causal_sample)<1.96*se_desc;
        in_desc_causal=abs(hat_theta-theta_causal)<1.96*se_desc;
        
        in_causal_sample_desc=abs(hat_theta-theta_desc)<1.96*se_causal_sample;
        in_causal_sample_causal_sample=abs(hat_theta-theta_causal_sample)<1.96*se_causal_sample;
        in_causal_sample_causal=abs(hat_theta-theta_causal)<1.96*se_causal_sample;
     
        in_causal_desc=abs(hat_theta-theta_desc)<1.96*se_causal;
        in_causal_causal_sample=abs(hat_theta-theta_causal_sample)<1.96*se_causal;
        in_causal_causal=abs(hat_theta-theta_causal)<1.96*se_causal;
          
        in_Z_causal_sample_desc=abs(hat_theta-theta_desc)<1.96*se_Z_causal_sample;
        in_Z_causal_sample_causal_sample=abs(hat_theta-theta_causal_sample)<1.96*se_Z_causal_sample;
        in_Z_causal_sample_causal=abs(hat_theta-theta_causal)<1.96*se_Z_causal_sample;
 
        in_Z_causal_desc=abs(hat_theta-theta_desc)<1.96*se_Z_causal;
        in_Z_causal_causal_sample=abs(hat_theta-theta_causal_sample)<1.96*se_Z_causal;
        in_Z_causal_causal=abs(hat_theta-theta_causal)<1.96*se_Z_causal;
 
        % Calculating sample coverage (ie, using estimated variances)
        in_hat_ehw_desc=abs(hat_theta-theta_desc)<1.96*se_hat_ehw;
        in_hat_ehw_causal_sample=abs(hat_theta-theta_causal_sample)<1.96*se_hat_ehw;
        in_hat_ehw_causal=abs(hat_theta-theta_causal)<1.96*se_hat_ehw;
 
        in_hat_desc_desc=abs(hat_theta-theta_desc)<1.96*se_hat_desc;
        in_hat_desc_causal_sample=abs(hat_theta-theta_causal_sample)<1.96*se_hat_desc;
        in_hat_desc_causal=abs(hat_theta-theta_causal)<1.96*se_hat_desc;
 
        in_hat_causal_sample_desc=abs(hat_theta-theta_desc)<1.96*se_hat_causal_sample;
        in_hat_causal_sample_causal_sample=abs(hat_theta-theta_causal_sample)<1.96*se_hat_causal_sample;
        in_hat_causal_sample_causal=abs(hat_theta-theta_causal)<1.96*se_hat_causal_sample;
 
        in_hat_causal_desc=abs(hat_theta-theta_desc)<1.96*se_hat_causal;
        in_hat_causal_causal_sample=abs(hat_theta-theta_causal_sample)<1.96*se_hat_causal;
        in_hat_causal_causal=abs(hat_theta-theta_causal)<1.96*se_hat_causal;
 
        in_hat_boot_desc=abs(hat_theta-theta_desc)<1.96*se_hat_boot;
        in_hat_boot_causal_sample=abs(hat_theta-theta_causal_sample)<1.96*se_hat_boot;
        in_hat_boot_causal=abs(hat_theta-theta_causal)<1.96*se_hat_boot;
  
        
        % Collecting output
        theta_out(isim,1:4)=[hat_theta-theta_desc,hat_theta-theta_causal_sample,hat_theta-theta_causal,hat_theta];
        
        out(isim,:)=[se_ehw,se_desc,se_causal_sample,se_causal,se_Z_causal_sample,se_Z_causal,se_hat_ehw,se_hat_desc,se_hat_causal_sample,se_hat_boot,se_hat_causal];
        
        in_ehw=[in_ehw_desc,in_ehw_causal_sample,in_ehw_causal];
        in_desc=[in_desc_desc,in_desc_causal_sample,in_desc_causal];
        in_causal_sample=[in_causal_sample_desc,in_causal_sample_causal_sample,in_causal_sample_causal];
        in_causal=[in_causal_desc,in_causal_causal_sample,in_causal_causal];
        in_Z_causal_sample=[in_Z_causal_sample_desc,in_Z_causal_sample_causal,in_Z_causal_sample_causal_sample];
        in_Z_causal=[in_Z_causal_desc,in_Z_causal_causal_sample,in_Z_causal_causal];
        in_hat_ehw=[in_hat_ehw_desc,in_hat_ehw_causal_sample,in_hat_ehw_causal];
        in_hat_desc=[in_hat_desc_desc,in_hat_desc_causal_sample,in_hat_desc_causal];
        in_hat_causal_sample=[in_hat_causal_sample_desc,in_hat_causal_sample_causal_sample,in_hat_causal_sample_causal];
        in_hat_causal=[in_hat_causal_desc,in_hat_causal_causal_sample,in_hat_causal_causal];
        in_hat_boot=[in_hat_boot_desc,in_hat_boot_causal_sample,in_hat_boot_causal];
        coverage(isim,:)=[in_ehw,in_desc,in_causal_sample,in_causal,in_Z_causal_sample,in_Z_causal,in_hat_ehw,in_hat_desc,in_hat_causal_sample,in_hat_boot,in_hat_causal];
        
        if floor(isim/csim)*csim==isim
           sd=[mean(theta_out(1:isim,1:3));std(theta_out(1:isim,1:3))]; 
           cov=mean(coverage(1:isim,:));
           mean_se=mean(out(1:isim,:));
           
           tabel(:,simulation_design)=[sd(2,:)';mean_se(1,1);cov(1,1:3)';mean_se(1,2);cov(1,4:6)';mean_se(1,3);cov(1,7:9)';mean_se(1,4);cov(1,10:12)';mean_se(1,7);cov(1,19:21)';mean_se(1,8);cov(1,22:24)';mean_se(1,9);cov(1,25:27)';mean_se(1,10);cov(1,28:30)';mean_se(1,11);cov(1,31:33)'];
           gtabel=gtable(tabel,3);
           save tabel_oct11 tabel gtabel simulation_design isim
           [simulation_design,isim];
        end
        
    end
    
 end
    
tabel;
gtabel=gtable(tabel,3);

% Saving Simulation Table
save tabel_oct11 tabel gtabel
