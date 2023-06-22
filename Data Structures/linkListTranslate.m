% =========================================================================
% FUNCTION linkListTranslate
% =========================================================================
function linkList = linkListTranslate(N, adjacency)

% Converts an adjacency matrix to a linkList by considering all
% possible neighbors a2 for agent a1. Adds a2 to the linkList for a1
% if the adjacency matrix specifies a connection.

linkList = [];

for a1 = 1:N
	linkList(a1,1) = 0;
	for a2 = 1:N
		neighborCheck = adjacency(a1, a2);
		if(neighborCheck == 1)
            % First index corresponds to number of neighbors
			linkList(a1,1) = 1 + linkList(a1,1);
			linkList(a1, linkList(a1,1) + 1) = a2;
		end
	end
end
end
