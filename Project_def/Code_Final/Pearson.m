% Pearson Correlation and CI
[corr_SG,~,RL_1,RU_1] = corrcoef(LGD,k_SG);
[corr_AR,~,RL_2,RU_2] = corrcoef(LGD,k_AR);

disp(['corr_SG = ',num2str(corr_SG(2,1))])
disp(['IC = (',num2str(RL_1(2,1)),', ',num2str(RU_1(2,1)),')'])
disp(['corr_AR = ',num2str(corr_AR(2,1))])
disp(['IC = (',num2str(RL_2(2,1)),', ',num2str(RU_2(2,1)),')'])

rho_Pearson = [corr_SG(2,1), corr_AR(2,1)];

% Scatter Plot
figure
scatter(k_AR,LGD)
lsline()
hold on
xlabel('k_{AR}')
ylabel('LGD')
title('All Rated Issuers')

figure
scatter(k_SG,LGD)
lsline()
hold on
xlabel('k_{SG}')
ylabel('LGD')
title('Speculative Grade Issuers')
