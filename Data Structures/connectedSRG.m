% =========================================================================
% FUNCTION connectedSRG
% =========================================================================

% =========================================================================
% The connectivity theorem used in this code is proved here:
% https://math.stackexchange.com/questions/864604/checking-connectivity-of-adjacency-matrix
%
% Alternate proof based on binomial theorem:
% https://math.stackexchange.com/questions/994847/adjacency-matrix-and-connectivity-proof
% =========================================================================

function [linkList, pos] = connectedSRG(N, r)

% Check for connectivity by noting that if connectivity(i, j) > 0 for each
% element, then the SRG is connected

connected = false;
while(~connected)
    [linkList, pos] = spatialRandomGraph(N,r);
    
    % Taking adjacency ^ N will represent the number of walks between 
    % agents i and j of length N while allowing self-loops; this provides
    % information on connectivity
    adjacency = adjacencyMatrix(N, linkList) + eye(N);
    
    connectivity = adjacency ^ N;
    if (all(connectivity(:) > 0))
        connected = true;
    end
end