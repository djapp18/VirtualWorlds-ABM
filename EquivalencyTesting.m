% =========================================================================
% Original Code Reference: ?Opinion polarization by learning from social
% feedback?,
% S. Banisch and E. Olbrich, The Journal of Mathematical Sociology, 2019,
% 43:2, 76-103
%
% =========================================================================

% =========================================================================
% Initializations
% =========================================================================

clear all; close all; clc;
tic;

% Default values for parameters taken from Banisch & Olbrich (2019)
N = 100;
r = 0.225;
max_steps = 20000 * N;

eps = 0.1;
alpha = 0.1;
beta = 1;
h = 0;

% =========================================================================
% Program configuration
% =========================================================================

% Determine whether to plot or write to file
single_run_plot_flag = 0;

if(single_run_plot_flag)
    write_file_flag = 0;
    iteration_max = 1;
    fprintf(['Single point plotting for visual equivalency testing ',...
        '.... \n\n']);    
else
    write_file_flag = 1;
    iteration_max = 100;
    fprintf(['Saving data to file for equivalency testing: Repeats = ',...
        '%d .... \n'], iteration_max);
end

% Open file to write out results for future processing
if (write_file_flag >0)
    outfile = 'cong_links_comparison_round_2.txt';
    file1 = fopen(outfile,'w');
    fprintf(file1, '%8s \n',[...
        'iteration      alpha      beta    r       h',...
        '      ratio_congruent_links    ratio_congruent_links_beta',...
        '      agents_with_opinion_difference']);
end

% Display parameters; can comment out to speed up run speed
if(single_run_plot_flag)
    % Print all the parameters to console
    fprintf(['[Agents N = %d ][r = %f][Steps = %d][eps =%f]',...
        '[alpha = %f][beta = %d][h = %f]\n'], N,r,max_steps, eps,alpha,beta,h);
end

% =========================================================================
% Repeat each data point iteration_max number of times
% =========================================================================

for iteration=1:1:iteration_max
    
    % Obtain SRG. Here, linkList represents the real world
    [linkList, pos] = connectedSRG(N,r);
    
    % Setup convictions, the first column corresponds to opinion "red" and the
    % second column corresponds to opinion "blue"
    Q = rand(N,2) - 0.5;
    Q_start = Q;
    Q_beta = Q;
    
    % =========================================================================
    % Simulation loop
    % =========================================================================
    
    for step = 1:max_steps
        % Get neighbor a2
        a1 = randi(N);
        std_rand = rand();
        
        numN = linkList(a1,1);
        neighborIndex = randi(numN);
        a2 = linkList(a1, neighborIndex+1);
        
        % =====================================================================
        % Model with eps
        % =====================================================================
        
        expression = (Q(a1,2) > Q(a1,1));
        if(std_rand <eps)
            expression = -(expression - 1);
        end
        
        % =====================================================================
        % Perform Q-learning update - with eps
        % =====================================================================
        
        reaction = (Q(a2,2) > Q(a2,1));
        reward = (expression * 2 - 1) * (reaction *2 - 1);
        
        Q(a1, expression + 1) = (1-alpha) * ...
            Q(a1, expression +1) + alpha * reward;
        
        % =====================================================================
        % Calculate metrics - with eps
        % =====================================================================
        
        % Dispersion and mean of convictions
        dQ = Q(:,1) - Q(:,2);
        sigS = var(dQ);
        dQBar = mean(dQ);
        
        sigSquare(step) = sigS;
        dQBar_step(step)= dQBar;
        
        % Congruent links
        ratio_congruent_links(step)= congruentLinks(N, linkList,Q);
        
        % =====================================================================
        % Model with beta
        % =====================================================================
        
        % Determine expression based on exploration rate. Here, an expression
        % of 0 corresponds to opinion "red" while an expression of 1 corresponds
        % to opinion "blue." The computation of eps is different here as the
        % numerator is hardcoded to be a single opinion; however, the relevant
        % expressions occur with the same probability as shown in the paper.
        eps_beta = exp( Q_beta(a1,2) * beta ) / ( exp( Q_beta(a1,1) * beta )...
            + exp( Q_beta(a1,2) * beta ) );
        if(std_rand() < eps_beta)
            expression_beta = 1;
        else
            expression_beta = 0;
        end
        
        % =====================================================================
        % Perform Q-learning update - with beta
        % =====================================================================
        
        reaction_beta = (Q_beta(a2,2) > Q_beta(a2,1));
        reward_beta = (expression_beta * 2 - 1) * (reaction_beta *2 - 1);
        
        Q_beta(a1, expression_beta + 1) = (1-alpha) * ...
            Q_beta(a1, expression_beta +1) + alpha * reward_beta;
        
        % =====================================================================
        % Calculate metrics - with beta
        % =====================================================================
        
        % Dispersion and mean of convictions
        dQ_beta = Q_beta(:,1) - Q_beta(:,2);
        sigS_beta = var(dQ_beta);
        dQBar_beta = mean(dQ_beta);
        
        sigSquare_beta(step) = sigS_beta;
        dQBar_step_beta(step)= dQBar_beta;
        
        % Congruent links
        ratio_congruent_links_beta(step)= congruentLinks(N, linkList,Q_beta);
        
    end
    
    % =========================================================================
    % Plotting and metrics
    % =========================================================================
    
    if (single_run_plot_flag)
        % Plot dispersion
        figure;plot((1:max_steps),sigSquare, (1:max_steps),sigSquare_beta);
        axis tight;grid on;ylim([0 4]);
        title('Mean dispersion over time');xlabel('Time steps');ylabel('Mean dispersion');
        
        % Plot congruent links
        figure;plot((1:max_steps),ratio_congruent_links, (1:max_steps),ratio_congruent_links_beta);
        axis tight;grid on;ylim([0.5 1]);
        title('Congruent links ratio over time');xlabel('Time steps');ylabel('Congruent links ratio');
        
        % Plot starting state, the final system state for both models
        drawSystemState(pos, linkList, Q_start);
        title('Starting State');
        
        drawSystemState(pos, linkList, Q);
        title('Constant \epsilon');
        
        drawSystemState(pos, linkList, Q_beta);
        title('Dynamic \beta');
    end
    
    % Determine similarity between two models by checking number of agents
    % with the same public opinion in both models
    agents_with_opn_difference = sum(abs((Q(:, 2) > Q(:, 1)) -...
        (Q_beta(:, 2) > Q_beta(:, 1))));
    
    fprintf(['[iteration = %d of %d ] [congruent_links = %f ]',...
        '[congruent_links_beta =  %f] [agents_with_opn_diff = %d]\n'],...
        iteration, iteration_max, ratio_congruent_links(step),...
        ratio_congruent_links_beta(step),agents_with_opn_difference);
    
    % Write relevant values to file
    if(write_file_flag > 0)
        fprintf(file1, ['%d          %8.3f   %8.3f %8.3f %8.3f    ',...
            '%8.6f                %8.6f                          %d\n'],...
            iteration, alpha, beta, r, h, ratio_congruent_links(step),...
            ratio_congruent_links_beta(step),agents_with_opn_difference);
    end
    
end

toc;
% =========================================================================
% Close file handle
% =========================================================================
if(write_file_flag > 0)
    fclose(file1);
end


