% This is a wrapper function to the MEX function stepFast.c

function stepMex(obj, nStep)

  % Select which step function to use
  stepFunc = [];
  
  if(obj.useLocalJumps)
    if(obj.typeFlag == 4)
      stepFunc = @obj.stepFastLocalChem;
    else
      stepFunc = @obj.stepFastLocal;
    end
  else
    % This function restricts selection of new synapses to neighbours
    % Note that stepFast (with obj.RGCwidth = 2, default value)
    % works similarly to stepFastGlobal. I.e. any RGC can connect
    % to a SC neuron, if RGCwidth less than 2 it becomes restricted
    
    % stepFunc = @obj.stepFast;
    stepFunc = @obj.stepFastGlobal;    
  end
      
  tic
    
  % Call the MEX function to step time nSteps

  try
    if(~exist('nStep'))
      nStep = obj.nSteps - obj.curStep;
    end

    idx = find(obj.curStep < obj.actionIteration ...
               & obj.actionIteration <= obj.curStep + nStep);
    
    
    if(~isempty(idx))
      
      shortSteps = obj.actionIteration(idx(1)) - obj.curStep;

      for i = 2:numel(idx)
        shortSteps(i) = obj.actionIteration(idx(i)) ...
                      - obj.actionIteration(idx(i-1));
      end
      
      shortSteps(end+1) = nStep - sum(shortSteps);
      
    else
      % No actions interferring - run all nStep at once
      shortSteps = nStep;
    end

    % Clear out the interrupts at distance 0
    shortSteps(find(shortSteps == 0)) = [];

    % Make sure there are no negative steps!!
    assert(all(shortSteps) > 0)
    
    for i = 1:numel(shortSteps)
  
      % We need to change from matlab to C indexing
      obj.presynapticConnections = obj.presynapticConnections - 1;

      for iM = 1:numel(obj.neighbourSC)
        obj.neighbourSC{iM} = obj.neighbourSC{iM} - 1;
      end

      for iM = 1:numel(obj.synapseNeighbourhood)
        obj.synapseNeighbourhood{iM} = obj.synapseNeighbourhood{iM} - 1;
      end
      % Done changing   
      
      [presynapticConnections, presynapticWeight, ...
       numPresynapticConnections, totalWeightRGC, ...
       totalWeightSC, time] = stepFunc(shortSteps(i)); 
      
      % We need to change back from C to matlab indexing
      for iM = 1:numel(obj.synapseNeighbourhood)
        obj.synapseNeighbourhood{iM} = obj.synapseNeighbourhood{iM} + 1;
      end
    
      for iM = 1:numel(obj.neighbourSC)
        obj.neighbourSC{iM} = obj.neighbourSC{iM} + 1;
      end

      obj.presynapticConnections = presynapticConnections + 1;
  
      obj.presynapticWeight = presynapticWeight;
      obj.numPresynapticConnections = numPresynapticConnections;
      obj.totalWeightRGC = totalWeightRGC;
      obj.totalWeightSC = totalWeightSC;
      
      obj.curStep = obj.curStep + shortSteps(i);
      obj.time = time;
      
      % Done converting back to matlab
      
      idx = find(obj.actionIteration == obj.curStep);
      
      for j = 1:numel(idx)
        fprintf('Performing action at iteration %d\n', ...
                obj.actionIteration(idx(j)))
        action = obj.actionHandler{idx(j)};
        
        if(~isempty(obj.actionParameters))
          parm = obj.actionParameters{idx(j)};
          action(parm);  
        else
          action()
        end
        
        disp('Action done.')
      end
      
    end
  catch e
    getReport(e)
    
    disp('Your data structures are currently in C-format')
    disp('Type return to resume stepMex, this will convert the indexes back')
    keyboard
    
  end
    
  
  
  unconnectedSC = nnz(obj.numPresynapticConnections == 0);
  unconnectedRGC = nnz(obj.totalWeightRGC == 0);
  
  if(unconnectedSC > 0)
    fprintf('%d SC with no presynaptic RGC\n', unconnectedSC)
  end
  
  if(unconnectedRGC > 0)
    fprintf('%d RGC with no postsynaptic SC\n', unconnectedRGC)
  end
  toc
  
  % Clear the post synaptic connection representation and the old
  % connection matrix, to make sure we do not mix old and new
  % connectivity
  
  obj.postsynapticConnections = [];
  obj.postsynapticWeight = []; 
  obj.connectionMatrix = [];
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
end