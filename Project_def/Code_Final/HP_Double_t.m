%% Naive approach
disp(['For the alpha: ', num2str(alpha*100), '%'])
disp(['For the nu: ',num2str(nu)])

[EL_SG_naive, RC_SG_naive]=Naive_Approach_tStudent(PD_SG_hat,LGD_hat,alpha,nu);
disp(['RC_SG_naive: ', num2str(RC_SG_naive)])

[EL_AR_naive, RC_AR_naive]=Naive_Approach_tStudent(PD_AR_hat,LGD_hat,alpha,nu);
disp(['RC_All_Rated_naive: ', num2str(RC_AR_naive)])
disp(' ')

%% Simulation
k=@(Select)norminv(Select);

factor=@(Select)(1-exp(-50*Select))/(1-exp(-50));
rho=@(Select)0.12*factor(Select)+0.24*(1-factor(Select));

rng(1)
N_sim=1e7;
I=50;
epsilon=trnd(nu,I,1);
M=trnd(nu,1,N_sim);
LGD_Simulated=std_LGD *trnd(nu,N_sim,1)+LGD_hat;
k_SG_Simulated=std_SG_k*trnd(nu,N_sim,1)+k_SG_hat;
k_AR_Simulated=std_AR_k*trnd(nu,N_sim,1)+k_AR_hat;

X=@(Select)sqrt(rho(Select))*M+sqrt(1-rho(Select))*epsilon;
X_SG=X(PD_SG_hat);
X_AR=X(PD_AR_hat);

%% Calling all the 8 (4+4) case of interest analysis
% k fix, LGD simulated, SG
[RC_k_SG,add_on_k_SG,EL_k_SG,VaR_k_SG]=add_on_Approach_tStudent_HP(0,LGD_Simulated,PD_SG_hat,M,X_SG,RC_SG_naive,EL_SG_naive,alpha,nu);
% k fix, LGD simulated, All Rated
[RC_k_AR,add_on_k_AR,EL_k_AR,VaR_k_AR]=add_on_Approach_tStudent_HP(0,LGD_Simulated,PD_AR_hat,M,X_AR,RC_AR_naive,EL_AR_naive,alpha,nu);
% k simulated, LGD fix, SG
[RC_LGD_SG,add_on_LGD_SG,EL_LGD_SG,VaR_LGD_SG]=add_on_Approach_tStudent_HP(1,k_SG_Simulated,LGD_hat,M,X_SG,RC_SG_naive,EL_SG_naive,alpha,nu);
% k simulated, LGD fix, All Rated
[RC_LGD_AR,add_on_LGD_AR,EL_LGD_AR,VaR_LGD_AR]=add_on_Approach_tStudent_HP(1,k_AR_Simulated,LGD_hat,M,X_AR,RC_AR_naive,EL_AR_naive,alpha,nu);
% k simulated, LGD simulated, SG
[RC_ind_SG,add_on_ind_SG,EL_ind_SG,VaR_ind_SG]=add_on_Approach_tStudent_HP(2,LGD_Simulated,k_SG_Simulated,M,X_SG,RC_SG_naive,EL_SG_naive,alpha,nu);
% k simulated, LGD simulated, All Rated
[RC_ind_AR,add_on_ind_AR,EL_ind_AR,VaR_ind_AR]=add_on_Approach_tStudent_HP(2,LGD_Simulated,k_AR_Simulated,M,X_AR,RC_AR_naive,EL_AR_naive,alpha,nu);
% k simulated, LGD simulated, SG Correlated
[LGD_Simulated_SG,k_SG_Simulated_SG]=Correlated_Distribution_tStudent(LGD_hat,std_LGD,k_SG_hat,std_SG_k,rho_Pearson(1),N_sim,nu);
[RC_corr_SG,add_on_corr_SG,EL_corr_SG,VaR_corr_SG]=add_on_Approach_tStudent_HP(2,LGD_Simulated_SG,k_SG_Simulated_SG,M,X_SG,RC_SG_naive,EL_SG_naive,alpha,nu);
% k simulated, LGD simulated, AR Correlated
[LGD_Simulated_AR,k_AR_Simulated_AR]=Correlated_Distribution_tStudent(LGD_hat,std_LGD,k_AR_hat,std_AR_k,rho_Pearson(2),N_sim,nu);
[RC_corr_AR,add_on_corr_AR,EL_corr_AR,VaR_corr_AR]=add_on_Approach_tStudent_HP(2,LGD_Simulated_AR,k_AR_Simulated_AR,M,X_AR,RC_AR_naive,EL_AR_naive,alpha,nu);
%% Display the result
AO=[add_on_k_AR,add_on_k_SG,add_on_LGD_AR,add_on_LGD_SG,add_on_ind_AR,add_on_ind_SG,add_on_corr_AR,add_on_corr_SG];
RC=[RC_k_AR,RC_k_SG,RC_LGD_AR,RC_LGD_SG,RC_ind_AR,RC_ind_SG,RC_corr_AR,RC_corr_SG];
VR=[VaR_k_AR,VaR_k_SG,VaR_LGD_AR,VaR_LGD_SG,VaR_ind_AR,VaR_ind_SG,VaR_corr_AR,VaR_corr_SG];
EL=[EL_k_AR,EL_k_SG,EL_LGD_AR,EL_LGD_SG,EL_ind_AR,EL_ind_SG,EL_corr_AR,EL_corr_SG];
Display_Result(AO,RC,VR,EL)
