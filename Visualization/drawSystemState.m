% =========================================================================
% FUNCTION drawSystemState
% =========================================================================
function[] = drawSystemState(pos,linkList,Q)
 
% Utilize agent positions, links, and opinions to generate plot revealing
% current state of the opinion dynamics system

figure;
set(gcf,'renderer','painters');
 
N = length(pos);
 
% Draw links as before        
for a1 = 1:N
    for a2 = 1:linkList(a1,1);
        line(   [pos(a1,1),pos(linkList(a1,a2+1),1)],...
                [pos(a1,2),pos(linkList(a1,a2+1),2)],...
                [1 1],'LineStyle','-','Color',[.7 .7 .7]);
    end;
    text(pos(a1,1),pos(a1,2),3,num2str(a1),'FontSize',12,...
            'HorizontalAlignment','center');
end;
hold on;

% Draw agents with s_i = 1
ncOut = (Q(:, 2) > Q(: ,1));
nodeSize = 1400*reshape(abs(Q(:, 2) - Q(: ,1)),N,1);
 
dQ = reshape(abs(Q(:, 2) - Q(: ,1)),N,1);
nsIn = (dQ/max(dQ)) .*  (0.66 * mean(nodeSize));   
        
ix = find((ncOut > 0));
scatter3(   pos(ix,1),pos(ix,2),zeros(length(ix),1)+2.2,1+nsIn(ix),...
            'MarkerEdgeColor','k','MarkerFaceColor','b',...
            'Marker','o','LineWidth',1);        
    
ix = find(ncOut == 0);
scatter3(   pos(ix,1),pos(ix,2),zeros(length(ix),1)+2.2,1+nsIn(ix),...
            'MarkerEdgeColor','k','MarkerFaceColor','r',...
            'Marker','o','LineWidth',1);
set(gcf,'Color','white')        
axis off;       
hold off;

end
