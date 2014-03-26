function calculateCU(obj)

  % assert(isempty(obj.SCz)); % Make sure SC is not 3D
  if(~isempty(obj.SCz))
    obj.calculateCU3D();
    return
  end
      
      
  disp('Calculating C and U, this can take some time.')
  
  % Lookup tables for activity contributions
  obj.CAct = exp(-obj.retinalDistance()/obj.bAct);
  
  % We are changing how the UAct is stored, instead of storing all
  % possible pairs, we only store them for each neighbour
  % obj.UAct = exp(-obj.SCDistance().^2/(2*obj.aAct^2));
  
  % OBS! stepFast assumes that UAct at distance 0 is 1
  
  obj.UAct = {};
  for i = 1:numel(obj.neighbourSC)
    obj.UAct{i} = zeros(size(obj.neighbourSC{i}));
    
    for j = 1:numel(obj.neighbourSC{i})
      idx = obj.neighbourSC{i}(j);
      d = sqrt((obj.SCap(i) - obj.SCap(idx))^2 ...
               + (obj.SCml(i) - obj.SCml(idx))^2);
      obj.UAct{i}(j) = exp(-d^2/(2*obj.aAct^2));
    end
    
  end
  

end