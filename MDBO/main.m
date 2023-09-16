clear all;
close all;
clc;

SearchAgents_no = 30;    % population size
Function_name = 'F1';    % Test Functions: F1-F23 
Max_iteration = 500;     % Maximum Iterations
cnt_max = 1;             % Runs

[lb, ub, dim, fobj] = Get_Functions_details(Function_name); % Load details of selected test function


Curve_DBO = zeros(1, Max_iteration);
Curve_MDBO = zeros(1, Max_iteration);

 
for cnt = 1 : cnt_max


    [DBO_Best_score(cnt),    DBO_Best_pos(cnt,:),   DBO_curve]   =   DBO(SearchAgents_no, Max_iteration, lb, ub, dim, fobj);      % DBO
    [MDBO_Best_score(cnt),   MDBO_Best_pos(cnt,:),  MDBO_curve]  =  MDBO(SearchAgents_no, Max_iteration, lb, ub, dim, fobj);      % MDBO
    
    
    Curve_DBO = Curve_DBO + DBO_curve;
	Curve_MDBO = Curve_MDBO + MDBO_curve;

end

% Iteration curve
Curve_DBO = Curve_DBO/cnt_max;
Curve_MDBO = Curve_MDBO/cnt_max;


% Std
std_DBO = std(DBO_Best_score);
std_MDBO = std(MDBO_Best_score);


% Worst
worst_DBO = max(DBO_Best_score);
worst_MDBO = max(MDBO_Best_score);


% Best
best_DBO = min(DBO_Best_score);
best_MDBO = min(MDBO_Best_score);


% Mean
mean_DBO = mean(DBO_Best_score);
mean_MDBO = mean(MDBO_Best_score);


figure(1)
semilogy(Curve_DBO,'color', [125, 195, 254]./255, 'linestyle', '-', 'linewidth', 2)
hold on
semilogy(Curve_MDBO,'color', [118, 80, 5]./255, 'linestyle', '-', 'linewidth', 1.5)

title(Function_name,'FontName', 'Helvetica', 'FontSize', 12);
xlabel('Iteration', 'FontName', 'Helvetica', 'FontSize', 12);
ylabel('Best score obtained so far', 'FontName', 'Helvetica', 'FontSize', 12);
set(gca,'linewidth',2);
grid on
box on
legend('DBO','MDBO')% Comparison of different metaheuristics
grid off


% Result Outputs
disp(['DBO£ºWorst: ', num2str(worst_DBO), ', Best: ', num2str(best_DBO), ', Mean: ', num2str(mean_DBO), ', Std: ', num2str(std_DBO)]);
disp(['MDBO£ºWorst: ', num2str(worst_MDBO), ', Best: ', num2str(best_MDBO), ', Mean: ', num2str(mean_MDBO), ', Std: ', num2str(std_MDBO)]);

