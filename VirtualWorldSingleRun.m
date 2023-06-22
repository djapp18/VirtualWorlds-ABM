% =========================================================================
% Original Code Reference: "Opinion polarization by learning from social
% feedback",
% S. Banisch and E. Olbrich, The Journal of Mathematical Sociology, 2019,
% 43:2, 76-103
%
% =========================================================================

% =========================================================================
% Initializations
% =========================================================================
  
% Random seed - REMOVE FOR TEST!
rng(1)

clear all; close all; clc;
tic;

% Default values for parameters taken from Banisch & Olbrich (2019)
N = 100;
r = 0.225;
max_steps = 20000 * N;
%max_steps = 200000;

alpha = 0.01;
beta = 1;
h = 4;

lambda = 0.4;
numWorlds = 2;

% Print all the parameters to console
fprintf('Agents N =  %d\n', N);
fprintf('r =  %f\n', r);
fprintf('Steps =  %d\n', max_steps);

fprintf('alpha = %f \n', alpha);
fprintf('beta = %d \n', beta);
fprintf('h =  %f\n', h);

fprintf('lambda =  %f \n', lambda);
fprintf('numWorlds =  %d \n', numWorlds);

% Obtain SRG. Here, linkList represents the real world
[linkList, pos] = connectedSRG(N,r);

% Setup convictions, the first column corresponds to opinion "red" and the
% second column corresponds to opinion "blue"
Q = rand(N,2) - 0.5;
Q_start = Q;

% Initialize virtual worlds.
theWorlds = virtualWorlds(N, numWorlds);

% =========================================================================
% Simulation loop
% =========================================================================

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
   
    % Default set of neighbors for a1 to interact with are contacts in the
    % real world
    numN = linkList(a1,1);
    linkListRow = linkList(a1, 2:(numN + 1));
    neighbors = linkListRow;
    
    % =========================================================================
    % Pick a virtual world
    % =========================================================================
    
    % Assume that there is add-1 smoothing for frequencies of both opinions
    % when calculating softmax; this prevents the possibility of
    % divide-by-zero errors
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

    % =========================================================================
    % Determine whether to enter virtual world
    % =========================================================================
    
    if (rand() < lambda)
        [theWorlds, virtualNeighbors] = theWorlds.login(a1, worldSelect);
        if(~isempty(virtualNeighbors))
            neighbors = virtualNeighbors;
            numN = length(neighbors);
        end
    end
    
    % =========================================================================
    % Incorporate homophily, pick a2 from possible neighbors
    % =========================================================================
    
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
    
    % =========================================================================
    % Perform Q-learning update
    % =========================================================================
    
    reaction = (Q(a2,2) > Q(a2,1));
    reward = (expression * 2 - 1) * (reaction *2 - 1);

    Q(a1, expression + 1) = (1-alpha) * ...
        Q(a1, expression +1) + alpha * reward;
    
    % =========================================================================
    % Calculate metrics
    % =========================================================================
    
    dQ = Q(:,1) - Q(:,2);
    
    % Dispersion and mean of convictions
    sigS = var(dQ);
    dQBar = mean(dQ);

    sigSquare(step) = sigS;
    dQBar_step(step)=dQBar;
    
    % Congruent links
    ratio_congruent_links(step)= congruentLinks(N, linkList,Q);
    
    % Virtual world metrics
    world = 1;
    worldOneParticipants = theWorlds.getParticipants(world);

    world = 2;
    worldTwoParticipants = theWorlds.getParticipants(world);

    % Proportion of opinions in subreddits
    numBlue = dot(worldOneParticipants, (Q(:,2) > Q(:,1))) + 1;
    numRed = dot(worldOneParticipants, (Q(:,2) <= Q(:,1))) + 1;
    worldOneRatios(step) = numBlue / numRed;

    numBlue = dot(worldTwoParticipants, (Q(:,2) > Q(:,1))) + 1;
    numRed = dot(worldTwoParticipants, (Q(:,2) <= Q(:,1))) + 1;
    worldTwoRatios(step) = numBlue / numRed;

    % Average convictions in each subreddit
    worldOneConvictions = Q(logical(worldOneParticipants), :);
    worldOnedQ = worldOneConvictions(:, 1) - worldOneConvictions(:, 2);
    worldOnedQBar = mean(worldOnedQ);
    worldOneAverages(step) = worldOnedQBar;

    worldTwoConvictions = Q(logical(worldTwoParticipants), :);
    worldTwodQ = worldTwoConvictions(:, 1) - worldTwoConvictions(:, 2);
    worldTwodQBar = mean(worldTwodQ);
    worldTwoAverages(step) = worldTwodQBar;
        
end

% =========================================================================
% Plotting and metrics
% =========================================================================

% Plot system states
drawSystemState(pos, linkList, Q_start);
title('Starting State');

drawSystemState(pos, linkList, Q);
title('Soft Max');

% Plot real world metrics
figure;plot((1:max_steps),sigSquare);
axis tight;grid on;ylim([0 4]);
title('Mean Dispersion');xlabel('time steps');ylabel('mean Dispersion');

figure;plot((1:max_steps),ratio_congruent_links);
axis tight;grid on;ylim([0.5 1]);
title('Congruent Links Ratio');xlabel('time steps');ylabel('ratio of congruent links');

% Plot virtual world metrics
figure;semilogy(1:max_steps,worldOneRatios,'-bo','LineWidth',2);
hold on;grid on;
semilogy(1:max_steps,worldTwoRatios,'-ro','LineWidth',2);
legend("Subreddit1", "Subreddit2");
title({['\lambda = ',num2str(lambda),'  \alpha = ',num2str(alpha)],...
    ['  h = ',num2str(h),'  steps = ',int2str(step),...
    '  r = ',num2str(r)],[]}, 'FontSize', 18);


% =========================================================================
% toc gives computational time
% =========================================================================
toc;
