%% Naive

[EL_SG_naive_V, RC_SG_naive_V] = Naive_Approach(PD_SG_hat,LGD_hat,alpha);
[EL_AR_naive_V, RC_AR_naive_V] = Naive_Approach(PD_AR_hat,LGD_hat,alpha);

for nu = 2:20
    [EL_SG_naive_t(nu-1), RC_SG_naive_t(nu-1)] = Naive_Approach_tStudent(PD_SG_hat,LGD_hat,alpha,nu);
    [EL_AR_naive_t(nu-1), RC_AR_naive_t(nu-1)] = Naive_Approach_tStudent(PD_AR_hat,LGD_hat,alpha,nu);
end

Naive_SG_V_99 = [RC_SG_naive_V; EL_SG_naive_V]
Naive_AR_V_99 = [RC_AR_naive_V; EL_AR_naive_V]
Naive_SG_t_99 = [RC_SG_naive_t; EL_SG_naive_t]
Naive_AR_t_99 = [RC_AR_naive_t; EL_AR_naive_t]

%% HP t-Student
for nu = 2:20
    factor_SG = (1-exp(-50*PD_SG_hat))/(1-exp(-50));
%     factor_AR = (1-exp(-50*PD_AR_hat))/(1-exp(-50));
    rho_SG = 0.12*factor_SG + 0.24*(1 - factor_SG );
%     rho_AR = 0.12*factor_AR + 0.24*(1 - factor_AR );
    k_SG = norminv(PD_SG_hat);
%     k_AR = norminv(PD_AR_hat);

    rng(1)
    N_sim = 2e4;
    I = 50;
    epsilon = trnd(nu,I,1);
    M = trnd(nu,1,N_sim);
    LGD_Simulated  = std_LGD *trnd(nu,N_sim,1) + LGD_hat;
    k_SG_Simulated = std_SG_k*trnd(nu,N_sim,1) + k_SG_hat;
%     k_AR_Simulated = std_AR_k*trnd(nu,N_sim,1) + k_AR_hat;

    X_SG = sqrt(rho_SG)*M + sqrt(1-rho_SG)*epsilon; 
%     X_AR = sqrt(rho_AR)*M + sqrt(1-rho_AR)*epsilon; 

    [EL_SG_naive, RC_SG_naive] = Naive_Approach_tStudent(PD_SG_hat,LGD_hat,alpha,nu);
%     [EL_AR_naive, RC_AR_naive] = Naive_Approach_tStudent(PD_AR_hat,LGD_hat,alpha,nu);

    [RC_k_SG(nu-1),add_on_k_SG(nu-1),EL_k_SG,VaR_k_SG] = add_on_Approach_tStudent_HP(0,LGD_Simulated,PD_SG_hat,M,X_SG,RC_SG_naive,EL_SG_naive,alpha,nu);
%     [RC_k_AR(nu-1),add_on_k_AR(nu-1),EL_k_AR,VaR_k_AR] = add_on_Approach_tStudent_HP(0,LGD_Simulated,PD_AR_hat,M,X_AR,RC_AR_naive,EL_AR_naive,alpha,nu);
    [RC_LGD_SG(nu-1),add_on_LGD_SG(nu-1),EL_LGD_SG,VaR_LGD_SG] = add_on_Approach_tStudent_HP(1,k_SG_Simulated,LGD_hat,M,X_SG,RC_SG_naive,EL_SG_naive,alpha,nu);
%     [RC_LGD_AR(nu-1),add_on_LGD_AR(nu-1),EL_LGD_AR,VaR_LGD_AR] = add_on_Approach_tStudent_HP(1,k_AR_Simulated,LGD_hat,M,X_AR,RC_AR_naive,EL_AR_naive,alpha,nu);
    [RC_ind_SG(nu-1),add_on_ind_SG(nu-1),EL_ind_SG,VaR_ind_SG] = add_on_Approach_tStudent_HP(2,LGD_Simulated,k_SG_Simulated,M,X_SG,RC_SG_naive,EL_SG_naive,alpha,nu);
%     [RC_ind_AR(nu-1),add_on_ind_AR(nu-1),EL_ind_AR,VaR_ind_AR] = add_on_Approach_tStudent_HP(2,LGD_Simulated,k_AR_Simulated,M,X_AR,RC_AR_naive,EL_AR_naive,alpha,nu);
   
    [LGD_Simulated_SG,k_SG_Simulated_SG] = Correlated_Distribution_tStudent(LGD_hat,std_LGD,k_SG_hat,std_SG_k,rho_Pearson(1),N_sim,nu);
    [RC_corr_SG(nu-1),add_on_corr_SG(nu-1),EL_corr_SG,VaR_corr_SG] = add_on_Approach_tStudent_HP(2,LGD_Simulated_SG,k_SG_Simulated_SG,M,X_SG,RC_SG_naive,EL_SG_naive,alpha,nu);
%     [LGD_Simulated_AR,k_AR_Simulated_AR] = Correlated_Distribution_tStudent(LGD_hat,std_LGD,k_AR_hat,std_AR_k,rho_Pearson(2),N_sim,nu);
%     [RC_corr_AR(nu-1),add_on_corr_AR(nu-1),EL_corr_AR,VaR_corr_AR] = add_on_Approach_tStudent_HP(2,LGD_Simulated_AR,k_AR_Simulated_AR,M,X_AR,RC_AR_naive,EL_AR_naive,alpha,nu);
end

k_SG_t_HP_99 = [RC_k_SG; add_on_k_SG]
LGD_SG_t_HP_99 = [RC_LGD_SG; add_on_LGD_SG]
ind_SG_t_HP_99 = [RC_ind_SG; add_on_ind_SG]
corr_SG_t_HP_99 = [RC_corr_SG; add_on_corr_SG]

k_AR_t_HP_99 = [RC_k_AR; add_on_k_AR]
LGD_AR_t_HP_99 = [RC_LGD_AR; add_on_LGD_AR]
ind_AR_t_HP_99 = [RC_ind_AR; add_on_ind_AR]
corr_AR_t_HP_99 = [RC_corr_AR; add_on_corr_AR]

%% LHP t-Student DA FAREEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
for nu = 2:20
%     factor_SG = (1-exp(-50*PD_SG_hat))/(1-exp(-50));
    factor_AR = (1-exp(-50*PD_AR_hat))/(1-exp(-50));
%     rho_SG = 0.12*factor_SG + 0.24*(1 - factor_SG );
    rho_AR = 0.12*factor_AR + 0.24*(1 - factor_AR );
%     k_SG = norminv(PD_SG_hat);
    k_AR = norminv(PD_AR_hat);

    rng(1)
    N_sim = 1e7;
    M = trnd(nu,N_sim,1);
    LGD_Simulated  = std_LGD *trnd(nu,N_sim,1) + LGD_hat;
%     k_SG_Simulated = std_SG_k*trnd(nu,N_sim,1) + k_SG_hat;
    k_AR_Simulated = std_AR_k*trnd(nu,N_sim,1) + k_AR_hat;

%     [EL_SG_naive, RC_SG_naive] = Naive_Approach_tStudent(PD_SG_hat,LGD_hat,alpha,nu);
    [EL_AR_naive, RC_AR_naive] = Naive_Approach_tStudent(PD_AR_hat,LGD_hat,alpha,nu);

%     [RC_k_SG(nu-1),add_on_k_SG(nu-1),EL_k_SG,VaR_k_SG] = add_on_Approach_tStudent(0,LGD_Simulated,PD_SG_hat,M,RC_SG_naive,EL_SG_naive,alpha,nu);
    [RC_k_AR(nu-1),add_on_k_AR(nu-1),EL_k_AR,VaR_k_AR] = add_on_Approach_tStudent(0,LGD_Simulated,PD_AR_hat,M,RC_AR_naive,EL_AR_naive,alpha,nu);
%     [RC_LGD_SG(nu-1),add_on_LGD_SG(nu-1),EL_LGD_SG,VaR_LGD_SG] = add_on_Approach_tStudent(1,k_SG_Simulated,LGD_hat,M,RC_SG_naive,EL_SG_naive,alpha,nu);
    [RC_LGD_AR(nu-1),add_on_LGD_AR(nu-1),EL_LGD_AR,VaR_LGD_AR] = add_on_Approach_tStudent(1,k_AR_Simulated,LGD_hat,M,RC_AR_naive,EL_AR_naive,alpha,nu);
%     [RC_ind_SG(nu-1),add_on_ind_SG(nu-1),EL_ind_SG,VaR_ind_SG] = add_on_Approach_tStudent(2,LGD_Simulated,k_SG_Simulated,M,RC_SG_naive,EL_SG_naive,alpha,nu);
    [RC_ind_AR(nu-1),add_on_ind_AR(nu-1),EL_ind_AR,VaR_ind_AR] = add_on_Approach_tStudent(2,LGD_Simulated,k_AR_Simulated,M,RC_AR_naive,EL_AR_naive,alpha,nu);
   
%     [LGD_Simulated_SG,k_SG_Simulated_SG] = Correlated_Distribution_tStudent(LGD_hat,std_LGD,k_SG_hat,std_SG_k,rho_Pearson(1),N_sim,nu);
%     [RC_corr_SG(nu-1),add_on_corr_SG(nu-1),EL_corr_SG,VaR_corr_SG] = add_on_Approach_tStudent(2,LGD_Simulated_SG,k_SG_Simulated_SG,M,RC_SG_naive,EL_SG_naive,alpha,nu);
    [LGD_Simulated_AR,k_AR_Simulated_AR] = Correlated_Distribution_tStudent(LGD_hat,std_LGD,k_AR_hat,std_AR_k,rho_Pearson(2),N_sim,nu);
    [RC_corr_AR(nu-1),add_on_corr_AR(nu-1),EL_corr_AR,VaR_corr_AR] = add_on_Approach_tStudent(2,LGD_Simulated_AR,k_AR_Simulated_AR,M,RC_AR_naive,EL_AR_naive,alpha,nu);
end

k_SG_t_LHP_99 = [RC_k_SG; add_on_k_SG]
LGD_SG_t_LHP_99 = [RC_LGD_SG; add_on_LGD_SG]
ind_SG_t_LHP_99 = [RC_ind_SG; add_on_ind_SG]
corr_SG_t_LHP_99 = [RC_corr_SG; add_on_corr_SG]

k_AR_t_LHP_99 = [RC_k_AR; add_on_k_AR]
LGD_AR_t_LHP_99 = [RC_LGD_AR; add_on_LGD_AR]
ind_AR_t_LHP_99 = [RC_ind_AR; add_on_ind_AR]
corr_AR_t_LHP_99 = [RC_corr_AR; add_on_corr_AR]

%% HP Vasicek
[EL_SG_naive, RC_SG_naive] = Naive_Approach(PD_SG_hat,LGD_hat,alpha);
[EL_AR_naive, RC_AR_naive] = Naive_Approach(PD_AR_hat,LGD_hat,alpha);

k = @(Select) norminv(Select);

factor = @(Select) (1-exp(-50*Select))/(1-exp(-50));
rho = @(Select) 0.12*factor(Select) + 0.24*(1 - factor(Select) );

rng(1)
N_sim = 2e4;
I = 50;
epsilon = randn(I,1);
LGD_Simulated  = std_LGD *randn(N_sim,1) + LGD_hat; % 20000x1
k_SG_Simulated = std_SG_k*randn(N_sim,1) + k_SG_hat;
k_AR_Simulated = std_AR_k*randn(N_sim,1) + k_AR_hat;
M = randn(1,N_sim);

X = @(Select) sqrt(rho(Select))*M + sqrt(1-rho(Select))*epsilon; % 50x20000
X_SG = X(PD_SG_hat);
X_AR = X(PD_AR_hat);

[RC_k_SG,add_on_k_SG,EL_k_SG,VaR_k_SG] = add_on_Approach_HP(0,LGD_Simulated,PD_SG_hat,M,X_SG,RC_SG_naive,EL_SG_naive,alpha);
[RC_k_AR,add_on_k_AR,EL_k_AR,VaR_k_AR] = add_on_Approach_HP(0,LGD_Simulated,PD_AR_hat,M,X_AR,RC_AR_naive,EL_AR_naive,alpha);
[RC_LGD_SG,add_on_LGD_SG,EL_LGD_SG,VaR_LGD_SG] = add_on_Approach_HP(1,k_SG_Simulated,LGD_hat,M,X_SG,RC_SG_naive,EL_SG_naive,alpha);
[RC_LGD_AR,add_on_LGD_AR,EL_LGD_AR,VaR_LGD_AR] = add_on_Approach_HP(1,k_AR_Simulated,LGD_hat,M,X_AR,RC_AR_naive,EL_AR_naive,alpha);
[RC_ind_SG,add_on_ind_SG,EL_ind_SG,VaR_ind_SG] = add_on_Approach_HP(2,LGD_Simulated,k_SG_Simulated,M,X_SG,RC_SG_naive,EL_SG_naive,alpha);
[RC_ind_AR,add_on_ind_AR,EL_ind_AR,VaR_ind_AR] = add_on_Approach_HP(2,LGD_Simulated,k_AR_Simulated,M,X_AR,RC_AR_naive,EL_AR_naive,alpha);

[LGD_Simulated_SG,k_SG_Simulated_SG] = Correlated_Distribution(LGD_hat,std_LGD,k_SG_hat,std_SG_k,rho_Pearson(1),N_sim);
[RC_corr_SG,add_on_corr_SG,EL_corr_SG,VaR_corr_SG] = add_on_Approach_HP(2,LGD_Simulated_SG,k_SG_Simulated_SG,M,X_SG,RC_SG_naive,EL_SG_naive,alpha);
[LGD_Simulated_AR,k_AR_Simulated_AR] = Correlated_Distribution(LGD_hat,std_LGD,k_AR_hat,std_AR_k,rho_Pearson(2),N_sim);
[RC_corr_AR,add_on_corr_AR,EL_corr_AR,VaR_corr_AR] = add_on_Approach_HP(2,LGD_Simulated_AR,k_AR_Simulated_AR,M,X_AR,RC_AR_naive,EL_AR_naive,alpha);

k_SG_V_HP_99 = [RC_k_SG; add_on_k_SG]
LGD_SG_V_HP_99 = [RC_LGD_SG; add_on_LGD_SG]
ind_SG_V_HP_99 = [RC_ind_SG; add_on_ind_SG]
corr_SG_V_HP_99 = [RC_corr_SG; add_on_corr_SG]

k_AR_V_HP_99 = [RC_k_AR; add_on_k_AR]
LGD_AR_V_HP_99 = [RC_LGD_AR; add_on_LGD_AR]
ind_AR_V_HP_99 = [RC_ind_AR; add_on_ind_AR]
corr_AR_V_HP_99 = [RC_corr_AR; add_on_corr_AR]

%% LHP Vasicek
rng(1)
N_sim = 10e7;
LGD_Simulated  = std_LGD *randn(N_sim,1) + LGD_hat;
k_SG_Simulated = std_SG_k*randn(N_sim,1) + k_SG_hat;
k_AR_Simulated = std_AR_k*randn(N_sim,1) + k_AR_hat;
M = randn(N_sim,1);

[EL_SG_naive, RC_SG_naive] = Naive_Approach(PD_SG_hat,LGD_hat,alpha);
[EL_AR_naive, RC_AR_naive] = Naive_Approach(PD_AR_hat,LGD_hat,alpha);

[RC_k_SG,add_on_k_SG,EL_k_SG,VaR_k_SG] = add_on_Approach(0,LGD_Simulated,PD_SG_hat,M,RC_SG_naive,EL_SG_naive,alpha);
[RC_k_SG_p,add_on_k_SG_p,EL_k_SG_p,VaR_k_SG_p] = add_on_Approach(0,mean(LGD_Simulated),PD_SG_hat,M,RC_SG_naive,EL_SG_naive,alpha);
[RC_k_AR,add_on_k_AR,EL_k_AR,VaR_k_AR] = add_on_Approach(0,LGD_Simulated,PD_AR_hat,M,RC_AR_naive,EL_AR_naive,alpha);
[RC_LGD_SG,add_on_LGD_SG,EL_LGD_SG,VaR_LGD_SG] = add_on_Approach(1,k_SG_Simulated,LGD_hat,M,RC_SG_naive,EL_SG_naive,alpha);
[RC_LGD_AR,add_on_LGD_AR,EL_LGD_AR,VaR_LGD_AR] = add_on_Approach(1,k_AR_Simulated,LGD_hat,M,RC_AR_naive,EL_AR_naive,alpha);
[RC_ind_SG,add_on_ind_SG,EL_ind_SG,VaR_ind_SG] = add_on_Approach(2,LGD_Simulated,k_SG_Simulated,M,RC_SG_naive,EL_SG_naive,alpha);
[RC_ind_AR,add_on_ind_AR,EL_ind_AR,VaR_ind_AR] = add_on_Approach(2,LGD_Simulated,k_AR_Simulated,M,RC_AR_naive,EL_AR_naive,alpha);

[LGD_Simulated_SG,k_SG_Simulated_SG] = Correlated_Distribution(LGD_hat,std_LGD,k_SG_hat,std_SG_k,rho_Pearson(1),N_sim);
[RC_corr_SG,add_on_corr_SG,EL_corr_SG,VaR_corr_SG] = add_on_Approach(2,LGD_Simulated_SG,k_SG_Simulated_SG,M,RC_SG_naive,EL_SG_naive,alpha);
[LGD_Simulated_AR,k_AR_Simulated_AR] = Correlated_Distribution(LGD_hat,std_LGD,k_AR_hat,std_AR_k,rho_Pearson(2),N_sim);
[RC_corr_AR,add_on_corr_AR,EL_corr_AR,VaR_corr_AR] = add_on_Approach(2,LGD_Simulated_AR,k_AR_Simulated_AR,M,RC_AR_naive,EL_AR_naive,alpha);

k_SG_V_LHP_99 = [RC_k_SG; add_on_k_SG]
LGD_SG_V_LHP_99 = [RC_LGD_SG; add_on_LGD_SG]
ind_SG_V_LHP_99 = [RC_ind_SG; add_on_ind_SG]
corr_SG_V_LHP_99 = [RC_corr_SG; add_on_corr_SG]

k_AR_V_LHP_99 = [RC_k_AR; add_on_k_AR]
LGD_AR_V_LHP_99 = [RC_LGD_AR; add_on_LGD_AR]
ind_AR_V_LHP_99 = [RC_ind_AR; add_on_ind_AR]
corr_AR_V_LHP_99 = [RC_corr_AR; add_on_corr_AR]
