%% Display of the basic statistical analysis on the dataset
disp(['min: ', num2str(round(min(LGD*100),2)), '| max: ', num2str(round(max(LGD*100),2)), '| mean: ',... 
    num2str(round(mean(LGD*100),2)), '| median: ', num2str(round(median(LGD*100),2)), '| std: ', ...
    num2str(round(std(LGD*100),2)), ' | obs: ', num2str(length(LGD))])

disp(['min: ', num2str(round(min(PD_SG*100),2)), ' | max: ', num2str(round(max(PD_SG*100),2)), '| mean: ', ... 
    num2str(round(mean(PD_SG*100),2)), '  | median: ', num2str(round(median(PD_SG*100),2)), ' | std: ', ...
    num2str(round(std(PD_SG*100),2)), '  | obs: ', num2str(length(PD_SG))])

disp(['min: ', num2str(round(min(PD_AR*100),2)), ' | max: ', num2str(round(max(PD_AR*100),2)), '    | mean: ',...
    num2str(round(mean(PD_AR*100),2)), ' | median: ', num2str(round(median(PD_AR*100),2)), ' | std: ', ...
    num2str(round(std(PD_AR*100),2)), '  | obs: ', num2str(length(PD_AR))])
disp(' ')

%% LGD distribution plot
figure
title('LGD distribution')
hold on
histfit(LGD,10)

%% Shapiro test
alpha = 0.95;

[H_LGD, pValue_LGD, W_LGD] = swtest(LGD, alpha);
[H_SG,  pValue_SG,  W_SG]  = swtest(k_SG, alpha);
[H_AR,  pValue_AR,  W_AR]  = swtest(k_AR, alpha);

disp(['LGD:  W: ', num2str(W_LGD), ' | p-value: ', num2str(pValue_LGD), ' | Normal distribution'])
disp(['k SG: W: ', num2str(W_SG), ' | p-value: ', num2str(pValue_SG), ' | Normal distribution'])
disp(['k AR: W: ', num2str(W_AR), ' | p-value: ', num2str(pValue_AR), ' | Normal distribution'])

%% Royston bivariate test
Roystest([LGD k_SG])
Roystest([LGD k_AR]) 

% qqplot graphs
figure
qqplot(LGD)
figure
qqplot(k_SG)
figure
qqplot(k_AR)
