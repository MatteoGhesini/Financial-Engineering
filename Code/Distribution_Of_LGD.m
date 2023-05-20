function [LGD_hat,std_LGD] = Distribution_Of_LGD(RR)
%[LGD_hat,std_LGD] = Distribution_Of_LGD(RR)
%The function computes the caracterizing parameter of the gaussian
%   distribution of LGD (Loss Given Default) starting from the observed RR
%   (Recovery Rates)
%
%INPUT
%   - RR: observed Recovery Rates
%
%OUTPUT
%   - LGD_hat: average of LGD distribution
%   - std_LGD: standard deviation of LGD distribution
%
%USING
%
%SUGGEST
%   
%

LGD = 1-RR;           % Computing Loss Given Default
LGD_hat = mean(LGD);    % Mean of LGD
std_LGD = std(LGD);     % Computing std

end