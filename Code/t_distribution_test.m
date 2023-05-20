function t_distribution_test(LGD,LGD_hat,k_SG,k_SG_hat,k_AR,k_AR_hat,nu)
%t_distribution_test(LGD,LGD_hat,k_SG,k_SG_hat,k_AR,k_AR_hat,nu)
%This function determines if two random samples, X1 and X2, are drawn from 
%the same underlying continuous population.
%
%INPUT
%   - LGD: vector of Loss Given Defaults
%   - LGD_hat: mean of LGD
%   - k_SG: vector of k Speculative Grade
%   - k_SG_hat: mean of k Speculative Grade
%   - k_AR: vector of k All Ratings
%   - k_AR_hat:mean of k All Ratings
%   - nu: Degrees of Freedom
%
%USING
%   - kstest2: to determine if two random samples are drawn from 
%     the same underlying continuous population.
%
%

LGD_test=kstest2(LGD, LGD_hat+std(LGD)*trnd(nu,length(LGD),1));     % Performing kstest2 between LGD and the simulated sample
k_SG_test=kstest2(k_SG, k_SG_hat+std(k_SG)*trnd(nu,length(k_SG),1));% Performing kstest2 between k_SG and the simulated sample
k_AR_test=kstest2(k_AR, k_AR_hat+std(k_AR)*trnd(nu,length(k_AR),1));% Performing kstest2 between k_AR and the simulated sample

% Display the results
if LGD_test==0
    disp(['LGD is a t-Student for nu: ',num2str(nu)])
else
    disp(['LGD is not a t-Student for nu: ',num2str(nu)])
end

if k_SG_test==0
    disp(['k_SG is a t-Student for nu: ',num2str(nu)])
else
    disp(['k_SG is not a t-Student for nu: ',num2str(nu)])
end

if k_AR_test==0
    disp(['k_AR is a t-Student for nu: ',num2str(nu)])
else
    disp(['k_AR is not a t-Student for nu: ',num2str(nu)])
end

end