function [RC,add_on,EL,VaR] = add_on_Approach(flag,simulated,fix_or_sim,M,RC_naive,EL_naive,alpha)
%[RC,add_on,EL,VaR] = add_on_Approach(flag,simulated,fix_or_sim,M,RC_naive,
% EL_naive,alpha)
%The function computes the new Regulatory Capital considering estimation on
%one(flag = 1,2) or two parameters (flag = 3), which are the Loss Given
%Default and the default point. The computation are in the case of Large
%Homogeneus Portfolio assumption and Vasiceck model, so the obligor is
%modeled as: X = sqrt(rho)*y + sqrt(1-rho)*e where y,e~N(0,1)
%
%INPUT
%   - flag: 0 k fix, LGD simulated
%           1 k simulated, LGD fix
%           2 k simulated, LGD simulated
%   - simulated: quantity that enter as vector previously simulated, for
%           each case of flag
%           0 LGD 
%           1 k 
%           2 LGD
%   - fix_or_sim: quantity that is fix (if we simulate only one parameter)
%           or simulated (otherwise) , for each case of flag
%           0 k 
%           1 LGD 
%           2 k
%   - M: vector simulated with market parameter
%   - RC_naive: Regulatory Capital in Naive approach 
%   - EL_naive: Expected Loss in Naive approach
%   - alpha: used in the computation of the VaR (usually 99.9%, or also
%           99%)
%
%OUTPUT
%   - RC: Regulatory Capital
%   - add_on: add-on
%   - EL: Expected Loss
%   - VaR: VaR at level alpha
%
%USING
%
%SUGGEST
%   - RC_naive, EL_naive are obtainable through Naive_approach function
%   - simulated, fix_or_sim (only for flag=2), M can be vector (column), in
%           that case they must have the same dimension
%
%

%% Error Section
if flag~=0 && flag~=1 && flag~=2
    disp('Flag not acceptable')
    return
end

if nargin<6
    disp('Too few input argument')
    return
end

if nargin<7
    disp('Using Default alpha 99.9%')
    alpha=0.999;
end
%%
TaylorApprox = @(k,PD) normcdf(k) - std(norminv(PD))^2/2 ... % Taylor approximation up to the third order  
    * k/sqrt(2*pi) * exp(-k^2/2) - PD;
switch flag
    case 0                                                   % Case: k fix, LGD simulated
        LGD=simulated;                                       % Simulation of LGD 
        PD=fix_or_sim;                                       % PD fixed
        k = fzero( @(k) TaylorApprox(k,PD) , 0);             % Computing k using fzero
        EL=mean(LGD)*PD;                                     % Computing EL
    case 1                                                   % Case: k simulated, LGD fix
        LGD=fix_or_sim;                                      % LGD fixed
        k=simulated;                                         % Simulation of k
        PD=normcdf(k);                                       % Computing PD 
        EL=LGD*mean(PD);                                     % Computing EL
    case 2                                                   % Case: k simulated, LGD simulated
        LGD=simulated;                                       % Simulation of LGD 
        k=fix_or_sim;                                        % k simulated 
        PD=normcdf(k);                                       % Computing PD 
        EL=mean(LGD)*mean(PD);                               % Computing EL
end

factor=(1-exp(-50*PD))/(1-exp(-50));                         
rho=0.12*factor+0.24*(1-factor);                             % rho computed as in Basel II 

EL_conditional=LGD.*normcdf((k-sqrt(rho).*M)./sqrt(1-rho));  % Computing EL conditional
VaR=prctile(EL_conditional,100*alpha);                       % VaR using prctile

RC=VaR-EL;                                                   % Computing Regulatory Capital
add_on=((RC-RC_naive)+(EL-EL_naive))/RC_naive;               % Obtaining the add_on

end