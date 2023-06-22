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
  
% Random seed - REMOVE FOR TEST!
% rng(1)

clear all; close all; clc;
tic;

% Default values for parameters taken from Banisch & Olbrich (2019)
N = 100;
r = 0.225;
max_steps = 20000 * N;
% max_steps = 2000000;

alpha = 0.01;
beta = 1;

numWorlds = 2;

% Print all the parameters to console
fprintf('Agents N =  %d\n', N);
fprintf('r =  %f\n', r);
fprintf('Steps =  %d\n', max_steps);

fprintf('alpha = %f \n', alpha);
fprintf('beta = %d \n', beta);

fprintf('numWorlds =  %d \n', numWorlds);

% Initialize virtual worlds.
theWorlds = virtualWorlds(N, numWorlds);

% Grid-search information
iteration_max = 100;

h_vector = [0 2 4 6];
lambda_vector = [0 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50];

write_file_flag = 1;

% =========================================================================
% Open file to write out results for future processing
% =========================================================================

if (write_file_flag >0)
    outfile = 'grid_search_lambda_h_NEW.txt';
    file1 = fopen(outfile,'w');
    fprintf(file1, '%11s \n', 'iteration      alpha      beta      lambda       h      consensus      Dispersion    dQbar    ratio_congruent_links    virtWorldDispRatio    virtWorldDispAverage');
end

% =========================================================================
% Three main loops: lambda - h grid search, where each data point is 
% repeated iteration_max number of times
% =========================================================================

for iteration=1:1:iteration_max
    % Repeat each data point iteration_max number of times
    fprintf('=== alpha = %f === beta = %d === iteration =  %d\n\n', alpha, beta, iteration);

    for lambda_index = 1:length(lambda_vector)
        % lambda loop - cycle through different values
        lambda = lambda_vector(lambda_index);
        fprintf('lambda =  %f\n', lambda);
        
        for h_index = 1:length(h_vector)
            % h loop - cycle through different values
            h = h_vector(h_index);
            fprintf('h =  %f\n', h);
            
            % Obtain SRG. Here, linkList represents the real world
            [linkList, pos] = connectedSRG(N,r);

            % Setup convictions, the first column corresponds to opinion 
            % "red" and the second column corresponds to opinion "blue"
            Q = rand(N,2) - 0.5;
            Q_start = Q;

            % Initialize virtual worlds.
            theWorlds = virtualWorlds(N, numWorlds);
            
            % =============================================================
            % Simulation loop
            % =============================================================
            for step = 1:max_steps

                a1 = randi(N);

                % Determine expression based on exploration rate. Here, an expression
                % of 0 corresponds to opinion "red" while an expression of 1 corresponds
                % to opinion "blue." The computation of eps is different here as the
                % numerator is hardcoded to be a single opinion; however, the relevant
                % expressions occur with the same probability as shown in the paper.
                eps = exp( Q(a1,2) * beta ) / ( exp( Q(a1,1) * beta ) + exp( Q(a1,2) * beta ) );
                if(rand() < eps)
                    expression = 1;
                else
                    expression = 0;
                end

                % Default set of neighbors for a1 to interact with are
                % contacts in the real world
                numN = linkList(a1,1);
                linkListRow = linkList(a1, 2:(numN + 1));
                neighbors = linkListRow;
                
                % =========================================================
                % Pick a virtual world
                % =========================================================
                
                % Assume that there is add-1 smoothing for frequencies of both 
                % opinions when calculating softmax; this prevents the possibility 
                % of divide-by-zero errors
                worldProbs = zeros(numWorlds, 1);
                for worldIndex = 1:numWorlds
                    theParticipants = theWorlds.getParticipants(worldIndex);
                    numBlue = dot(theParticipants, (Q(:,2) > Q(:,1))) + 1;
                    numRed = dot(theParticipants, (Q(:,2) <= Q(:,1))) + 1;

                    % Determination of the group based on intrinsic conviction (i.e.,
                    % Q(a1, 2) > Q(a1, 1))
                    if (Q(a1, 2) > Q(a1, 1))
                        ratio = numBlue/numRed;
                    else
                        ratio = numRed/numBlue;
                    end

                    worldProbs(worldIndex, :) = exp(ratio);
                end

                worldSelect = randsample(2, 1, true, worldProbs);
                
                % =========================================================
                % Determine whether to enter virtual world
                % =========================================================

                if (rand() < lambda)
                    [theWorlds, virtualNeighbors] = theWorlds.login(a1, worldSelect);
                    if(~isempty(virtualNeighbors))
                        neighbors = virtualNeighbors;
                        numN = length(neighbors);
                    end
                end

                % =========================================================
                % Incorporate homophily, pick a2 from possible neighbors
                % =========================================================

                % Expression for homophily taken from Maes & Bischofberger (2015)
                probs = zeros(numN, 1);
                dq1 = Q(a1,2) - Q(a1,1);

                for agentIndex = 1:numN
                    agent = neighbors(agentIndex);
                    dq2 = Q(agent, 2) - Q(agent, 1);

                    % Extra factor of 1/2 required, dq1 and dq2 range from -2 to 2
                    similarity = 1 - (abs(dq1 - dq2) / 4);        
                    factor = exp(similarity * h);        
                    probs(agentIndex) = factor;
                end

                neighborIndex = randsample(numN, 1, true, probs);
                a2 = neighbors(neighborIndex);
                
                % =========================================================
                % Perform Q-learning update
                % =========================================================

                reaction = (Q(a2,2) > Q(a2,1));
                reward = (expression * 2 - 1) * (reaction *2 - 1);

                Q(a1, expression + 1) = (1-alpha) * ...
                    Q(a1, expression +1) + alpha * reward;

            end
            
            % =============================================================
            % Calculate metrics
            % =============================================================
            
            dQ = Q(:,1) - Q(:,2);

            % Dispersion and mean of convictions
            sigS = var(dQ);
            dQBar = mean(dQ);

            % Congruent links
            ratio_congruent_links = congruentLinks(N, linkList,Q);
            
            % Consensus metrics
            
            consensus = false;
            consensusCheck = Q(:, 1) > Q(:, 2);
            
            szOpinionOne = sum(consensusCheck);
            szOpinionTwo = N - szOpinionOne;
            
            if (szOpinionOne > szOpinionTwo)
                szLargest = szOpinionOne;
            else
                szLargest = szOpinionTwo;
            end
            
            if(szOpinionOne == N || szOpinionTwo == N)
                consensus = true;
            end

            % Virtual world metrics
            world = 1;
            worldOneParticipants = theWorlds.getParticipants(world);
            
            world = 2;
            worldTwoParticipants = theWorlds.getParticipants(world);
            
            % Proportion of opinions in subreddits
            numBlue = dot(worldOneParticipants, (Q(:,2) > Q(:,1))) + 1;
            numRed = dot(worldOneParticipants, (Q(:,2) <= Q(:,1))) + 1;
            worldOneRatio = numBlue / numRed;

            numBlue = dot(worldTwoParticipants, (Q(:,2) > Q(:,1))) + 1;
            numRed = dot(worldTwoParticipants, (Q(:,2) <= Q(:,1))) + 1;
            worldTwoRatio = numBlue / numRed;
            
            virtWorldDispRatio = var([log10(worldOneRatio) log10(worldTwoRatio)]);

            % Average convictions in each subreddit
            worldOneConvictions = Q(logical(worldOneParticipants), :);
            worldOnedQ = worldOneConvictions(:, 1) - worldOneConvictions(:, 2);
            worldOnedQBar = mean(worldOnedQ);

            worldTwoConvictions = Q(logical(worldTwoParticipants), :);
            worldTwodQ = worldTwoConvictions(:, 1) - worldTwoConvictions(:, 2);
            worldTwodQBar = mean(worldTwodQ);

            % Edge case: if one virtual world is empty, default to zero
            % dispersion
            if(~sum(worldOneParticipants) || ~sum(worldTwoParticipants))
                virtWorldDispAverage = 0;
            end
            
            virtWorldDispAverage = var([worldOnedQBar worldTwodQBar]);
            
            % Save dispersion and mean values for all the iterations to
            % plot at the end of the run for visualization purpose
            super_dispersion(iteration,lambda_index,h_index) = ...
                sigS;
            super_dQBar(iteration,lambda_index,h_index) = dQBar;
            
            super_conLinksRatio(iteration,lambda_index,h_index) = ratio_congruent_links;
            super_virtWorldDispRatio(iteration,lambda_index,h_index) = virtWorldDispRatio;
            super_virtWorldDispAverage(iteration,lambda_index,h_index) = virtWorldDispAverage;
            
            if(consensus)
                super_consensus(iteration,lambda_index,h_index) = ...
                    1;
            else
                super_consensus(iteration,lambda_index,h_index) = ...
                    0;
            end
            
            % Write relevant values to file
            if(write_file_flag > 0)
                fprintf(file1, '%d           %8.3f  %8.3f    %8.3f    %8.3f    %8.3f     %8.3f      %8.3f     %8.3f      %8.3f        %8.3f \n',...
                    iteration, alpha, beta, lambda, h, consensus*1, sigS, dQBar,...
                    ratio_congruent_links, virtWorldDispRatio, virtWorldDispAverage);
            end
            
        end
        
    end
end

% =========================================================================
% Close file handle
% =========================================================================
if(write_file_flag > 0)
    fclose(file1);
end

% =========================================================================
% toc gives computational time
% =========================================================================
toc;
