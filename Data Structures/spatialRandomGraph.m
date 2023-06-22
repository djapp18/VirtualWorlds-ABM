% =========================================================================
% FUNCTION spatialRandomGraph
% =========================================================================
function[linkList, pos] = spatialRandomGraph(N,r)

% Generate all positions randomly, then for each agent a1 determine if 
% every other agent a2 is within r distance. If so, then add a2 to the 
% linkList for a1.

pos = rand(N,2);
linkList = [];

for a1 = 1:N
	linkList(a1,1) = 0;
	for a2 = 1:N
		dist = sqrt((pos(a1,1) - pos(a2,1))^2+ (pos(a1,2) - pos(a2,2))^2);
		if(ne(a1,a2) && dist < r)
            % First index corresponds to number of neighbors
			linkList(a1,1) = 1 + linkList(a1,1);
			linkList(a1, linkList(a1,1) + 1) = a2;
		end
	end
end
end