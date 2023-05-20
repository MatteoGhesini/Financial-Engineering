function [LGD_Simulated,k_Simulated] = Correlated_Distribution_tStudent(LGD_hat,std_LGD,k_hat,std_k,rho,N_sim,nu)
%[LGD_Simulated,k_Simulated] = Correlated_Distribution_tStudent(LGD_hat,
% std_LGD,k_hat,std_k,rho,N_sim,nu)
%The function computes the t-Student correlated distribution given two
%t-Student distribution and their correlation
%
%INPUT
%   - mu1: mean of the first t-Student distribution
%   - std1: standard deviation of the first t-Student distribution 
%   - mu2: mean of the second t-Student distribution
%   - std2: standard deviation of the second t-Student distribution
%   - rho: correlation between the two distribution
%   - N_sim: how many value to simulate of the correlated distribution
%
%OUTPUT
%   - Simulated1: first variable simulated considering correlation 
%   - Simulated2: second variable simulated considering correlation 
%
%USING
%   - my_mvtrand
%
%SUGGEST
%   - As rho it could be used the Pearson correlation
%   - Check before if the two marginal are t-Student through a Two-sample 
%       Kolmogorov-Smirnov goodness-of-fit hypothesis test between your
%       data and a simulated t-Student (with correct mean and standard
%       deviation)
%   - For the internal simulation we use rng(1), so we suggest to use it
%       also in your code
%
%

rng(1)                                  % Setting the seed
mu=[LGD_hat k_hat];                     % Vector of means
Cov=rho*std_LGD*std_k;                  % Computing the covariance 
Sigma=[std_LGD^2 Cov; Cov std_k^2];     % Computing variance covariance matrix

sample=my_mvtrand(N_sim,mu',Sigma,nu);  % Obtaining random vectors from the multivariate normal distribution

LGD_Simulated=sample(:,1);              % Taking the first column of sample
k_Simulated=sample(:,2);                % Taking the second column of sample

end