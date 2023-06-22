% =========================================================================
% FUNCTION congruentLinks
% =========================================================================
function congruentLinksRatio = congruentLinks(N,linkList,Q)
% This function finds the ratio of congruent links to total links 

cLinks = 0;
totalLinks = 0;

for a1 = 1:N
    length = linkList(a1, 1);
    for index = 1:length
        a2 = linkList(a1, index + 1);

        expression = (Q(a1,2) > Q(a1,1));
        reaction = (Q(a2,2) > Q(a2,1));

        % Only increment cLinks if expression and reaction are congruent
        if (expression == reaction)
            cLinks = cLinks + 1;
        end

        totalLinks = totalLinks + 1;
    end
end

cLinks = cLinks / 2;
totalLinks = totalLinks / 2;
congruentLinksRatio = cLinks/totalLinks;