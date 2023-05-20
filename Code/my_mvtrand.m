function distribution = my_mvtrand(N_sim,mu,sigma,nu,seed)
%distribution = my_mvtrand(N_sim,mu,sigma,nu)
%This function returns random numbers of Multivariate t distribution
% 
%INPUT
%   - N_sim: number of random values you want to simulate
%   - mu: pX1 vector of bivariate t distribution
%   - sigma: pXp Sigma matrix
%   - nu: Degree of Freedom
%   - seed: random seed for the simulations (default = 1)
% 
%OUTPUT
%   - distribution: nXp random numbers drawn from the bivariate t 
%       distribution
%
%USING
%   - gamrnd: for the generation of Gamma values
%
%

if nargin<5
    seed=1;                         % Setting the default seed 
end

rng(seed);
p = length(mu);                     % Initialization of p

Y=gamrnd(nu/2,1,[N_sim 1])/(nu/2);  % Simulation of a Chi-square
G=1./Y;                             % Computing Inverse-Chi-square
Z=mvnrnd(zeros(1,p),sigma,N_sim);   % Obtaining random vectors from the multivariate normal distribution

distribution=mu'+sqrt(G).*Z;        % Obtaining random numbers drawn from the bivariate t distribution

end