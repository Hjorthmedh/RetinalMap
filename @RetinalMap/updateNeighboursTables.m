function updateNeighboursTables(obj)

  disp('Updating neighbour tables')
  
  % RGC axons close to respective SC
  
  if(obj.useLocalJumps)
    if(isnan(obj.maxConnections))
      disp('You must set maxConnections')
      keyboard
    end
  else
    % We are doing global jumps
      
    % We only need space for this many connections
    % (the connections can have a weight > 1, so can have more synapses)
    if(isnan(obj.maxConnections))
      obj.maxConnections = numel(obj.nRGC);
    else
      obj.maxConnections = min(obj.nRGC,obj.maxConnections);
    end
  
  end
  
  fprintf('Using maxConnections = %d\n', obj.maxConnections)
  
  % SC neurons within interaction distance
  % exp(-r^2/(2*a^2)), at 5*a distance 3.7e-6 remains
  % at 4*a distance 3.5e-4 remains
  
  obj.neighbourSC = {};
  obj.nNeighbourSC = zeros(obj.nSC,1,'int32');
  dMax = obj.aAct*4;
  
  fprintf('Cutting off SC interactions after %f times decay width\n', ...
          dMax/obj.aAct);
  
  % These are all other SC neurons within "interaction range"
  
  for i = 1:obj.nSC
    
    if(isempty(obj.SCz))
      % SC is 2D
      idx = find(sqrt((obj.SCap(i)-obj.SCap).^2 ...
                      +(obj.SCml(i)-obj.SCml).^2) <= dMax);
    else
      % SC is 3D
      idx = find(sqrt((obj.SCap(i)-obj.SCap).^2 ...
                      +(obj.SCml(i)-obj.SCml).^2 ...
                      +(obj.SCz(i)-obj.SCz).^2) ...                      
                 <= dMax);   
    end
    
    % Update: Exclude the neuron self
    idx = setdiff(idx,i);

    obj.neighbourSC{i} = int32(idx);
    obj.nNeighbourSC(i) = numel(idx);
    
  end
  
  
  % SC neurons that the axon can jump to if it already synapses to
  % the i:th neuron
  
  obj.synapseNeighbourhood = {};
  obj.nSynapseNeighbourhood = zeros(obj.nSC,1,'int32');
  
  if(obj.useLocalJumps)
  
    fprintf('For local computation, using max synapse jump length %.4f\n',...
            obj.maxSynapseJumpLength)
  
    for i = 1:obj.nSC
      if(isempty(obj.SCz))
        % 2D SC
        idx = find(sqrt((obj.SCap(i)-obj.SCap).^2 ...
                        +(obj.SCml(i)-obj.SCml).^2) <= obj.maxSynapseJumpLength);
      else
        % 3D SC
        idx = find(sqrt((obj.SCap(i)-obj.SCap).^2 ...
                        +(obj.SCml(i)-obj.SCml).^2 ...
                        +(obj.SCz(i)-obj.SCz).^2 ) ...
                   <= obj.maxSynapseJumpLength);
        
      end
      
      % %  Update: Exclude the neuron self
      % idx = setdiff(idx,i);

      obj.synapseNeighbourhood{i} = int32(idx);
      obj.nSynapseNeighbourhood(i) = numel(idx);
  
    end
  end
  
end
