function importConnectionMatrix(obj,connectionFile,stateFile)

  if(~exist('stateFile') | isempty(stateFile))
    disp('You should load the old simulation state file to get gradients,')
    disp('and neuron positions for retina and SC.')
  else
    obj.loadState(stateFile);
  end  
  
  fprintf('Loading connection matrix from (assuming %d rows (SC), %d columns (RGC))\n%s\n', ...
          obj.nSC, obj.nRGC, connectionFile)
  try
    % Catherines format
    obj.connectionMatrix = dlmread(connectionFile,' ', ...
                                   [1 1 obj.nRGC obj.nSC]);
    obj.connectionMatrix = transpose(obj.connectionMatrix);
  catch
    % Fallback method
    % Stephens format for Gierer 2D
    obj.connectionMatrix = load(connectionFile);
  end
  

  
  obj.convertConnectionTables('mat2pre');
  obj.convertConnectionTables('pre2post');
  
  
  
  obj.totalWeightRGC = zeros(obj.nRGC,1);
  obj.totalWeightSC = zeros(obj.nSC,1);
  
  for i = 1:obj.nRGC
    nCon = obj.numPostsynapticConnections(i);
    obj.totalWeightRGC(i) = sum(obj.postsynapticWeight(1:nCon,i));
  end
  
  for i = 1:obj.nSC
    nCon = obj.numPresynapticConnections(i);
    obj.totalWeightSC(i) = sum(obj.presynapticWeight(1:nCon,i));    
  end
end