function [RC,add_on,EL,VaR] = add_on_Approach_HP(flag,simulated,fix_or_sim,M,X,RC_naive,EL_naive,alpha,GraphFlag)
%[RC,add_on,EL,VaR] = add_on_Approach_HP(flag,simulated,fix_or_sim,M,X,
% RC_naive,EL_naive,alpha)
%The function computes the new Regulatory Capital considering estimation on
%one(flag = 1,2) or two parameters (flag = 3), which are the Loss Given
%Default and the default point. The computation are in the case of
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
%   - X: matrix simulated with obligor states as by 
%           X = sqrt(rho)*M + sqrt(1-rho)*e
%           with M,e ~ N(0,1)
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
%The input matrix X is built as 
%   X = sqrt(rho)*M + sqrt(1-rho)*e    with M,e ~ N(0,1)
%   M: 1xN_simulation
%   e: Ix1 -> X: IxN_simulation
%   Doing so, X will be a matrix combining the I (number of obligor) case
%   with the N_simulation describing N different market condition (common
%   to all obligor)
%
%SUGGEST
%   - RC_naive, EL_naive are obtainable through Naive_approach function
%   - simulated, fix_or_sim (only for flag=2), M can be vector
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

if nargin<9
    GraphFlag=0;
end
%%
TaylorApprox = @(k,PD) normcdf(k) - std(norminv(PD))^2/2 ... % Taylor approximation up to the third order
    * k/sqrt(2*pi) * exp(-k^2/2) - PD;
switch flag
    case 0                                                   % Case: k fix, LGD simulated
        LGD=mean(simulated);                                 % Mean of simulated LGD
        PD=fix_or_sim;                                       % PD fixed
        k = fzero( @(k) TaylorApprox(k,PD) , 0);             % Computing k using fzero
    case 1                                                   % Case: k simulated, LGD fix
        LGD=fix_or_sim;                                      % LGD fixed
        k=mean(simulated);                                   % Mean of simulated k
        PD=normcdf(k);                                       % Computing PD
    case 2                                                   % Case: k simulated, LGD simulated
        LGD=mean(simulated);                                 % Mean of simulated LGD
        k=mean(fix_or_sim);                                  % Mean of k
        PD=normcdf(k);                                       % Computing PD
end

factor=(1-exp(-50*PD))/(1-exp(-50)); 
rho=0.12*factor+0.24*(1-factor );                            % rho computed as in Basel II 

%%
defaults=(X<k);                                              % Matrix with 1 in (i,j) if the jth obligor 
                                                             % will default in the ith market situation
defaulted_at_time=sum(defaults)';                            % Vector containing how many obligors 
                                                             % defaulted in a certain market condition
I=size(X,1);                                                 % Number of obligors considered
P_y=@(y)normcdf((k-sqrt(rho).*y)/sqrt(1-rho));               
if I>50
    my_nchoosek=@(I,m)1./(sqrt(2*pi*m.*(1-m./I)).*...        % Binomial distribution
        (m./I).^m.*(1-m./I).^(I-m));
else
    my_nchoosek=@(I,m)nchoosek(I,m);
end
P_m_y=@(m,y)my_nchoosek(I,m)*(P_y(y).^m).*(1-P_y(y)).^(I-m);                    

%% Loss Distribution 
% Vector presenting in each market condition how would be our loss,
% considering how many obligor defaulted, at step all the obligor have the
% same LGD
Loss=sum((X<k).*LGD')';

yGraph=zeros(I,1);                                           % Used to make a plot
new_Loss=zeros(I,1);                                         % Empirical Loss of m defaults
EL=zeros(I,1);                                               % Empirical Expected Loss of m defaults
for ii=1:I
    if ~isnan(find(defaulted_at_time==ii)) 
        f=find(defaulted_at_time==ii);                       % Search how many times (and when) there are a certain number of defaults
        new_Loss(ii)=Loss(f(1));                             % Save only one time the value of the loss for m defaults
        yGraph(ii)=size(f,1)/size(X,2);
        EL(ii)=new_Loss(ii)*yGraph(ii);                      % Loss*PD where PD is computed as how many times we see that number of defaults over all the simulation
    else
        new_Loss(ii)=0;                                      % Otherwise it would put NaN (case in which we never see that number of defaults)
        EL(ii)=0;                                            % Otherwise it would put NaN (case in which we never see that number of defaults)
        yGraph(ii)=0;                                        % Used to make a plot
    end
end
% now we have a vector of length equal to I, each element of which is the
% loss of having that number of default (e.g. new_Loss(3) = Loss|3defaults)
EL=mean(EL);

% "Analitical" expected loss of m defaults conditional
EL_conditional=zeros(I,1);
for ii=1:I
    EL_conditional(ii)=new_Loss(ii)*mean(P_m_y(ii,M),2);     % Computing EL conditional
end
%%
VaR=prctile(EL_conditional,100*alpha);                       % VaR using prctile
RC=VaR-EL;                                                   % Computing Regulatory Capital
add_on=((RC-RC_naive)+(EL-EL_naive))./RC_naive;              % Obtaining the add_on

%% To have a good plot we suggest to use 500 obligors (we have already 
% implemented che Kullback-Leibler approximation to speed up the 
% computation in that case)
if GraphFlag == 1
    VaR_999=prctile(EL_conditional,100*0.999);
    VaR_99=prctile(EL_conditional,100*0.99);
    figure
    plot(smoothdata(new_Loss/I),smoothdata(yGraph),'b')
    hold on
    xline(EL,'r-',{'EL'},'LabelOrientation','horizontal','LabelVerticalAlignment','middle','LineWidth',1)
    xline(VaR_99,'r-',{'VaR_{99%}'},'LabelVerticalAlignment','middle','LabelHorizontalAlignment','left','LineWidth',1)
    xline(VaR_999,'r-',{'VaR_{99.9%}'},'LabelVerticalAlignment','middle','LineWidth',1)
    title('Loss HP')
end
end