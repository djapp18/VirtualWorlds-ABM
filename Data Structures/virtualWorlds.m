% =========================================================================
% CLASS virtualWorlds
% =========================================================================
classdef virtualWorlds

% Each virtual world is similar to the real world except the first column
% now represents whether the specified agent is in the virtual world or not

properties (Access = private)
    worlds;
    numAgents;
end

methods
    % =====================================================================
    %  Constructor
    % =====================================================================
    function obj = virtualWorlds(N, numWorlds)
        obj.worlds = zeros(N, N + 1, numWorlds);
        obj.numAgents = N;
    end
    
    % =====================================================================
    % Returns a vector detailing which agents are in the virtual world
    % corresponding to worldIndex; if no index is provided, a matrix
    % corresponding to all of the virtual worlds is returned
    % =====================================================================
    function participants = getParticipants(obj, worldIndex)     
        if nargin < 1
            participants = squeeze(obj.worlds(:, 1, :));
        else
            participants = squeeze(obj.worlds(:, 1, worldIndex));
        end
    end
    
    % =====================================================================
    % Returns a matrix detailing the system state corresponding to
    % worldIndex. The info is organized the same as the real world. If no
    % index is provided, a 3D matrix corresponding to all of the virtual
    % worlds is returned.
    % =====================================================================
    function systemState = getSystemState(obj, worldIndex)     
        if nargin < 1
            systemState = squeeze(obj.worlds(:, 2:end, :));
        else
            systemState = squeeze(obj.worlds(:, 2:end, worldIndex));
        end
    end
    
    % =====================================================================
    % Agent a1 attempts to login into the virtual world corresponding to
    % worldSelect. If they are the first user to login, returns null.
    % =====================================================================
    function [obj, neighbors] = login(obj, a1, worldSelect)
        world = obj.worlds(:, :, worldSelect);
        
        % If agent a1 has not visited this virtual world before, they need
        % to "sign up"
        if (world(a1, 1) == 0)
            world(a1, 1) = 1;

            % Add this agent as a neighbor to all other agents in the
            % virtual world
            row = 1:obj.numAgents;
            row(row == a1) = [];
            world(row, 2) = world(row, 2) + 1;
            col = world(row, 2) + 2;
            idx = sub2ind(size(world), row, transpose(col));
            world(idx) = a1;

        end

        % If agent does not have any neighbors in virtual world (i.e., they
        % are the first to sign up) then they default back to their real
        % world neighbors.
        if(world(a1, 2) ~= 0)
            numN = world(a1, 2);
            neighbors = world(a1, 3:(numN + 2));
        else
            neighbors = [];
        end

        % Update general virtualWorlds data structure
        obj.worlds(:, :, worldSelect) = world; 
    end
end

end