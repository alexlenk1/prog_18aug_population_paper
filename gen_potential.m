
function [Z,xi,theta]=gen_population(n,psi,gamma,lambda,sig_theta,sig_xi,KZ);

% generating population data

%input
%n = population size
%psi = parameter determining the dependence between Z and theta
%gamma = parameter determining the dependence between Z and xi
%lambda = parameter determining the dependence between Z and U
%sig_theta = sd of theta 
%sig_xi = sd of xi
%KZ = dimension of Z, excluding constant term

%output 
%Z = n by (KZ+1) matrix of fixed attributes for population n
%xi = the intercept vector of size n in the potential outcome function
%theta = treatment effect vector of size n  in the potential outcome function

Z=[ones(n,1),randn(n,KZ)];  

xi_error=randn(n,1)*sig_xi; % xi_error is N(0, (sig_xi)^2 )
theta_error=randn(n,1)*sig_theta; % theta_error is N(0, (sig_theta)^2 )

xi=Z*gamma+xi_error;  % intercept for potential outcomes, distributed as N(Z*gamma, var(xi_error) ) 
theta=Z*psi+theta_error;      % treatment effect, distributed as N(Z*psi, var(theta_error) ) 
