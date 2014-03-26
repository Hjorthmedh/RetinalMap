function importConnectionIndexes(obj,indexFile,weightFile,stateFile)

  if(~exist('stateFile') | isempty(stateFile))
    disp('You should load the old simulation state file to get gradients,')
    disp('and neuron positions for retina and SC.')
  else
    obj.loadState(stateFile);
  end  
  
  % Data files has one row per RGC, multiple columns with indexes
  % to SC neurons (or weights).
  
  % We want one column per RGC, so transpose 
  index = transpose(load(indexFile));
  
  if(~exist('weightFile')  | isempty(weightFile))
    disp('No weight file given, assuming all weights are 1.');
    weight = ones(size(index));
  else
    weight = transpose(load(weightFile));
  end
  
  obj.postsynapticConnections = index;
  obj.postsynapticWeight = weight;
  obj.numPostsynapticConnections = sum(weight > 0,1);
  obj.totalWeightRGC = transpose(sum(weight,1));
  
  % Just some sanity checks on the input
  assert(size(obj.postsynapticConnections,2) == obj.nRGC);
  assert(max(obj.postsynapticConnections(:)) <= obj.nSC);
  assert(all(obj.postsynapticWeight(:)) >= 0);
  
  obj.convertConnectionTables('post2pre');
  
  obj.totalWeightSC = transpose(sum(obj.presynapticWeight,1));
  
end