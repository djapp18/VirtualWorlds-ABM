% =========================================================================
% FUNCTION adjacencyMatrix
% =========================================================================
function adjacency = adjacencyMatrix(N, linkList)

% To make adjacency matrix from linkList, create N x N of zeros, then scan
% through linkList for given agent a1. For each neighbor a2 specified, 
% replace 0 with 1 in adjacency matrix (~O(N^2) time)

adjacency = zeros(N);

for a1 = 1:N
    length = linkList(a1, 1);
    for index = 1:length
        a2 = linkList(a1, index + 1);
        adjacency(a1, a2) = 1;
    end
end
end