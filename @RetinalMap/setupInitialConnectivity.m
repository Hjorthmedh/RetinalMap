function setupInitialConnectivity(obj)

 
  % Pre allocate
  obj.numPresynapticConnections = zeros(obj.nSC,1); % Keep it double, for rand
  obj.totalWeightSC = zeros(obj.nSC,1);
  obj.totalWeightRGC = zeros(obj.nRGC,1);

  % Place the retinal axons
  % These have no impact unless synapse placement is restricted to
  % near existing synapses, in which case there need to be initial
  % synapses placed (with or without bias)
  if(obj.useLocalJumps)
    
    if(obj.RGCmlSpread == inf | isnan(obj.RGCmlSpread))
      obj.randomizeRGCAxonPositionNoBias();
    else
      obj.randomizeRGCAxonPositionWithBias();
    end   
    
  end
  
  % This lists all possible neighbours, and counts them
  obj.updateNeighboursTables();
  
  % Allocate space so we can connect to all neighbours if needed
  obj.presynapticConnections = ...
      zeros(obj.maxConnections,obj.nSC,'int32');
  obj.presynapticWeight = zeros(obj.maxConnections,obj.nSC);

  obj.calculateCU();
  
  % Place initial synapses (only for local jumps)
  obj.placeInitialSynapses(); 

  
end
