%% Run_GROUP_1A
% Ghesini Matteo, Toschi Anna

clear all
close all
clc

%% Dataset uploading
if ispc()                                                                      % Windows version
    [~, PD_SG, PD_AR, RR] = readExcel('Data\CreditModelRisk_RawData.xlsx');
else                                                                           % MacOS version - for semplicity we load the data previously saved
    load('Data\PD_AR.mat')                                                     % Loading data AR
    load('Data\PD_SG.mat')                                                     % Loading data SG
    load('Data\RR.mat')                                                        % Loading data RR
end

%% Normality Assumptions
LGD = 1-RR;                                                                    % Loss Given Default
k_SG = norminv(PD_SG);                                                         % Computing k_SG
k_AR = norminv(PD_AR);                                                         % Computing k_AR

%% Main statistical analysis
Statistica                                                                     % Statistical analysis
Pearson                                                                        % Computing Pearson correlation and its 95% confidence interval

%% Distributions initialization
[LGD_hat,std_LGD] = Distribution_Of_LGD(RR);                                   % Computing the mean and std of LGD
disp(['LGD ~ N(',num2str(LGD_hat*100),'%, ',num2str(std_LGD*100),'^2%)'])

[PD_SG_hat, k_SG_hat, std_SG_k] = Distribution_Of_k(PD_SG);                    % Computing the mean and std of k_SG
disp(['k_SG ~ N(',num2str(k_SG_hat),', ',num2str(std_SG_k),'^2)'])

[PD_AR_hat, k_AR_hat, std_AR_k] = Distribution_Of_k(PD_AR);                    % Computing the mean and std of k_AR
disp(['k_All_Rated ~ N(',num2str(k_AR_hat),', ',num2str(std_AR_k),'^2)'])
disp(' ')

%% Select the alpha
% alpha = 0.99;
alpha = 0.999;

%%
% LHP Approach 
LHP_Vasicek                                                                    % Computing LHP approach

% HP Approach 
HP_Vasicek                                                                     % Computing  HP approach

%% 
nu = 20;                                                                       % Varying nu parameter to have the different cases
t_distribution_test(LGD,LGD_hat,k_SG,k_SG_hat,k_AR,k_AR_hat,nu)                % Check if they are t-Student using ks-test2

% LHP Approach 
LHP_Double_t                                                                   % Computing LHP approach via Double t-student

% HP Approach
HP_Double_t                                                                    % Computing HP approach via Double t-student

% end
