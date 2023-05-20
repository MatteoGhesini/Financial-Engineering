function [PD_hat, k_hat, std_k] = Distribution_Of_k(PD)
%[PD_hat, k_hat, std_k] = Distribution_Of_k(PD)
%The function computes the caracterizing parameter of the gaussian
%   distribution of k (default point) starting from the observed PD
%
%INPUT
%   - PD: observed Probabilities to Default
%
%OUTPUT
%   - PD_hat: average of PD
%   - k_hat: average of k distribution
%   - std_k: standard deviation of k distribution
%
%USING
%
%SUGGEST
%   - k_hat it's not a simple mean but it's obtained as solution of
%  PD_hat = normcdf(k_hat) - std_k^2/2 * k_hat/sqrt(2*pi) * exp(-k_hat^2/2) 
%
%

PD_hat = mean(PD);                                  % Computing the mean of PD

k = norminv(PD);                                    % Computing k
std_k = std(k);                                     % Computing std 

myfun = @(k_hat) normcdf(k_hat) - std_k^2/2 *...
    k_hat/sqrt(2*pi) * exp(-k_hat^2/2) - PD_hat;
k_hat = fzero( @(k_hat) myfun(k_hat) , 0);          % Computing k_hat using fzero

end