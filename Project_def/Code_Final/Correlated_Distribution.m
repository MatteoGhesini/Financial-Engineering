function [Simulated1,Simulated2] = Correlated_Distribution(mu1,std1,mu2,std2,rho,N_sim)
%[Simulated1,Simulated2] = Correlated_Distribution(mu1,std1,mu2,std2,rho,
% N_sim)
%The function computes the normal correlated distribution given two normal
%distribution and their correlation
%
%INPUT
%   - mu1: mean of the first gaussian distribution
%   - std1: standard deviation of the first gaussian distribution 
%   - mu2: mean of the second gaussian distribution
%   - std2: standard deviation of the second gaussian distribution
%   - rho: correlation between the two distribution
%   - N_sim: how many value to simulate of the correlated distribution
%
%OUTPUT
%   - Simulated1: first variable simulated considering correlation 
%   - Simulated2: second variable simulated considering correlation 
%
%USING
%   - mvnrnd
%
%SUGGEST
%   - As rho it could be used the Pearson correlation
%   - Check before if the two marginal are gaussian through a Shapiro-Wilks
%       test
%   - Check before if the bivariate distribution is gaussian through a
%       Royston test (extended Shapiro-Wilks)
%   - For the internal simulation we use rng(1), so we suggest to use it
%       also in your code
%
%

rng(1)                                        % Setting the seed
mu=[mu1 mu2];                                 % Vector of means
Cov=rho*std1*std2;                            % Computing the covariance 
Sigma=[std1^2 Cov; Cov std2^2];               % Computing variance covariance matrix

sample=mvnrnd(mu,Sigma,N_sim);                % Obtaining random vectors from the multivariate normal distribution

Simulated1=sample(:,1);                       % Taking the first column of sample
Simulated2=sample(:,2);                       % Taking the second column of sample

end