function [EL_naive, RC_naive] = Naive_Approach_tStudent(PD_hat,LGD_hat,alpha,nu)
%[EL_naive, RC_naive] = Naive_Approach_tStudent(PD_hat,LGD_hat,alpha,nu)
%The function computes the Regulatory Capital considering the Naive IRB 
%Approach. The computation are in the case of Double t-Student model, so 
%the obligor is modeled as: X =sqrt(rho)*y + sqrt(1-rho)*e where y,e~t(0,1)
%
%INPUT
%   - PD_hat: Scalar mean of PD that we observe in the past
%   - LGD_hat: Scalar mean of LGD that we observe in the past 
%   - alpha: used in the computation of the VaR (usually 99.9%, or also
%           99%)
%   - nu: Degree of Freedom for the t-Student
%
%OUTPUT
%   - EL: Expected Loss
%   - RC: Regulatory Capital
%
%USING
%
%SUGGEST
%
%

factor=(1-exp(-50*PD_hat))/(1-exp(-50)); 
rho=0.12*factor+0.24*(1-factor );                       % Computing rho as in Basel II

EL_naive=LGD_hat*PD_hat;                                % Computing EL with Naive approach
k=norminv(PD_hat);                                      % Computing k 
RC_naive=LGD_hat*tcdf((k-sqrt(rho)*tinv(1-alpha,nu))... % Computing RC with Naive approach via
                /sqrt(1-rho),nu)-EL_naive;              % Double t-Student model

end