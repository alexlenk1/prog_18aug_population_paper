% main program, April 10th, 2014
clear   % ensure that the workspace is empty

load seed_oct2

application=false;    % false if you want to run simulations,
                      % true if you want to run the application with the
                      % nls data
                      
Nboot=1000;           % number of bootstrap replications

%  application to nls data
if application
   load nls_2012.txt
   luwe=nls_2012(:,1);
   educ=nls_2012(:,2);
   exper=nls_2012(:,3);
   age=nls_2012(:,4);
   fed=nls_2012(:,5);
   moed=nls_2012(:,6);
   kww=nls_2012(:,7);
   iq=nls_2012(:,8);
   clear nls_2012;
   
   U=educ;
   Z=[ones(length(educ),1),age,fed,moed,kww,iq];
   Y=luwe;
   beta=inv([U Z]'*[U Z])*([U Z]'*Y)
   theta=beta(1,1);
   
   se_ehw=se_ehw_calc(Y,U,Z);
   se_z=se_z_calc(Y,U,Z);
   [lcb_boot,ucb_boot,se_boot]=se_boot_calc(Y,U,Z,Nboot)
   % point estimate, ehw standard errors, proposed se_z standard errors,
   % bootstrap
   'est se_ehw se_z se_boot'
   [theta,se_ehw,se_z,se_boot]

   [mean(luwe),size(luwe)]
end      % end of nls application

% parameters for simulations
Nsim=50000;      % number of simulations. Should correspond to number of simulations in table
csim=1000;       % number of simulations before output is sent to screen for monitoring. Not important.
EN=1000;         % expected number of observations in the sample.
                 % if the population size is n, and the sampling rate is
                 % rho, EN is approximately rho*n. EN and rho are fixed in the
                 % simulations, and the population size is calculated as
                 % EN/rho.
KZ=2;            % number of fixed characteristics Z, beyond intercept.

tabel=zeros(39,3);   % matrix containing output for tables
                     % there are seven simulation designs

for simulation_design=1:7,
    if simulation_design==1,
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
    
    if simulation_design==2,
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

    if simulation_design==3,
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
    
    if simulation_design==4,
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
    
    if simulation_design==5,
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
    
    if simulation_design==6,
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
    
    if simulation_design==7,
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
    % the potential outcome function is
    % y(x)=xi_i+theta_i*x
    [Z,xi,theta]=gen_potential(n,psi,gamma,lambda,sig_theta,sig_xi,KZ);
    
    % calculate theta_causal and Omega
    theta_causal=mean(phi); % correct this
    Omega=zeros(KZ+3,KZ+3);
    
    theta_out=zeros(Nsim,4);
    out=zeros(Nsim,11);
    coverage=zeros(Nsim,33);
    
    [Z,alpha,phi];
    
    Delta_ehw=1+3*(psi'*psi+sig_eps*sig_eps);
    Delta_Z=Delta_ehw-(psi'*psi);
    Delta_cond=Delta_ehw-(psi'*psi+sig_eps*sig_eps);
    %pause
    H=1;               % variance of X, so with X randn(n,1), this is equal to 1
    V_ehw=inv(H)*Delta_ehw*inv(H)/(rho*n);
    V_desc=(1-rho)*inv(H)*Delta_ehw*inv(H)/(rho*n);
    V_causal_sample=inv(H)*Delta_cond*inv(H)/(rho*n);
    V_causal=inv(H)*(rho*Delta_cond+(1-rho)*Delta_ehw)*inv(H)/(rho*n);
    V_Z_causal_sample=inv(H)*Delta_Z*inv(H)/(rho*n);
    V_Z_causal=inv(H)*(rho*Delta_Z+(1-rho)*Delta_ehw)*inv(H)/(rho*n);
    
    se_ehw=sqrt(V_ehw);
    se_desc=sqrt(V_desc);
    se_causal_sample=sqrt(V_causal_sample);
    se_causal=sqrt(V_causal);
    se_Z_causal_sample=sqrt(V_Z_causal_sample);
    se_Z_causal=sqrt(V_Z_causal);
    
    for isim=1:Nsim
        [YR,XR,ZR,W,UR,R,TOmega,theta_desc,theta_causal_sample]=gen_sample(alpha,phi,Z,rho,xi);
        
        theta_causal_sample=R'*phi/sum(R);
        
        beta=inv([UR ZR]'*[UR ZR])*([UR ZR]'*YR);
        beta=inv([XR ZR]'*[XR ZR])*([XR ZR]'*YR);
        hat_theta=beta(1,1);    % calculate theta estimator
        
        
        
        se_hat_ehw=se_ehw_calc(YR,UR,ZR);
        se_hat_desc=sqrt(1-rho)*se_hat_ehw;
        se_hat_Z=se_z_calc(YR,UR,ZR);
        se_hat_causal_Z=se_causal_z_calc(YR,UR,ZR,rho);
        [lcb_boot,ucb_boot,se_hat_boot]=se_boot_calc(YR,UR,ZR,Nboot);
        
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
 
        in_hat_ehw_desc=abs(hat_theta-theta_desc)<1.96*se_hat_ehw;
        in_hat_ehw_causal_sample=abs(hat_theta-theta_causal_sample)<1.96*se_hat_ehw;
        in_hat_ehw_causal=abs(hat_theta-theta_causal)<1.96*se_hat_ehw;
 
        in_hat_desc_desc=abs(hat_theta-theta_desc)<1.96*se_hat_desc;
        in_hat_desc_causal_sample=abs(hat_theta-theta_causal_sample)<1.96*se_hat_desc;
        in_hat_desc_causal=abs(hat_theta-theta_causal)<1.96*se_hat_desc;
 
        in_hat_Z_desc=abs(hat_theta-theta_desc)<1.96*se_hat_Z;
        in_hat_Z_causal_sample=abs(hat_theta-theta_causal_sample)<1.96*se_hat_Z;
        in_hat_Z_causal=abs(hat_theta-theta_causal)<1.96*se_hat_Z;
 
        in_hat_causal_Z_desc=abs(hat_theta-theta_desc)<1.96*se_hat_causal_Z;
        in_hat_causal_Z_causal_sample=abs(hat_theta-theta_causal_sample)<1.96*se_hat_causal_Z;
        in_hat_causal_Z_causal=abs(hat_theta-theta_causal)<1.96*se_hat_causal_Z;
 
        in_hat_boot_desc=abs(hat_theta-theta_desc)<1.96*se_hat_boot;
        in_hat_boot_causal_sample=abs(hat_theta-theta_causal_sample)<1.96*se_hat_boot;
        in_hat_boot_causal=abs(hat_theta-theta_causal)<1.96*se_hat_boot;
  
        
        % collecting output
        theta_out(isim,1:4)=[hat_theta-theta_desc,hat_theta-theta_causal_sample,hat_theta-theta_causal,hat_theta];
        
        out(isim,:)=[se_ehw,se_desc,se_causal_sample,se_causal,se_Z_causal_sample,se_Z_causal,se_hat_ehw,se_hat_desc,se_hat_Z,se_hat_boot,se_hat_causal_Z];
        
        in_ehw=[in_ehw_desc,in_ehw_causal_sample,in_ehw_causal];
        in_desc=[in_desc_desc,in_desc_causal_sample,in_desc_causal];
        in_causal_sample=[in_causal_sample_desc,in_causal_sample_causal_sample,in_causal_sample_causal];
        in_causal=[in_causal_desc,in_causal_causal_sample,in_causal_causal];
        in_Z_causal_sample=[in_Z_causal_sample_desc,in_Z_causal_sample_causal,in_Z_causal_sample_causal_sample];
        in_Z_causal=[in_Z_causal_desc,in_Z_causal_causal_sample,in_Z_causal_causal];
        in_hat_ehw=[in_hat_ehw_desc,in_hat_ehw_causal_sample,in_hat_ehw_causal];
        in_hat_desc=[in_hat_desc_desc,in_hat_desc_causal_sample,in_hat_desc_causal];
        in_hat_Z=[in_hat_Z_desc,in_hat_Z_causal_sample,in_hat_Z_causal];
        in_hat_causal_Z=[in_hat_causal_Z_desc,in_hat_causal_Z_causal_sample,in_hat_causal_Z_causal];
        in_hat_boot=[in_hat_boot_desc,in_hat_boot_causal_sample,in_hat_boot_causal];
        coverage(isim,:)=[in_ehw,in_desc,in_causal_sample,in_causal,in_Z_causal_sample,in_Z_causal,in_hat_ehw,in_hat_desc,in_hat_Z,in_hat_boot,in_hat_causal_Z];
        
        if floor(isim/csim)*csim==isim,
           sd=[mean(theta_out(1:isim,1:3));std(theta_out(1:isim,1:3))] 
           cov=mean(coverage(1:isim,:))
           mean_se=mean(out(1:isim,:))
           
           tabel(:,simulation_design)=[sd(2,:)';mean_se(1,1);cov(1,1:3)';mean_se(1,2);cov(1,4:6)';mean_se(1,3);cov(1,7:9)';mean_se(1,4);cov(1,10:12)';mean_se(1,7);cov(1,19:21)';mean_se(1,8);cov(1,22:24)';mean_se(1,9);cov(1,25:27)';mean_se(1,10);cov(1,28:30)';mean_se(1,11);cov(1,31:33)']
           gtabel=gtable(tabel,3)
           save tabel_oct11 tabel gtabel simulation_design isim
           [simulation_design,isim]
        end
        
    end,
    
 end
    
tabel
gtabel=gtable(tabel,3)


save tabel_oct11 tabel gtabel
