function loadState(obj, fileName, skipTables)

  if(~exist('fileName') | isempty(fileName))
    fileName = sprintf('%s/%s.mat', obj.dataPath, obj.simName);
  end    
    
  fprintf('Loading state from %s\n', fileName)
  save = load(fileName);

  allFields = fieldnames(obj); 
  excludeFields = {'UAct','CAct','neighbourRGC','neighbourSC'};
  loadFields = setdiff(allFields,excludeFields);

  for i = 1:numel(loadFields)
    try
      eval(sprintf('obj.%s = save.old.%s;', ...
                   loadFields{i}, ...
                   loadFields{i}));
    catch 
      fprintf('Unable to load field %s\n', loadFields{i})
    end
  end
  
  % Recalculate the neighbour tables, this saves us space relative
  % to storing them on disk.
  
  if(~exist('skipTables'))
    skipTables = 1;
  end
  
  switch(skipTables)
    case 0
      disp('Recalculating neighbour relations')
      obj.updateNeighboursTables();
  
      % Recalculate C and U
      obj.calculateCU();
      
    case 1
      disp('Skipping neighbour tables, C and U.')

    case 2
      disp('Memory conservation: Clearing all but connection data')
      
      obj.RGCEphA = [];
      obj.RGCEphB = [];
      obj.RGCephrinA = [];
      obj.RGCephrinB = [];
      
      obj.RGCml = [];
      
      obj.SCephrinA = [];
      obj.SCephrinB = [];
      obj.SCEphA = [];
      obj.SCEphB = [];
      
      obj.SCmask = [];
      
      obj.nNeighbourSC = [];
      obj.nNeighbourRGC = [];
      obj.nSynapseNeighbourhood = [];
      
      obj.maxConnections = max(obj.numPresynapticConnections);
      obj.presynapticConnections = obj.presynapticConnections(1:obj.maxConnections,:);
      
      obj.presynapticWeight = obj.presynapticWeight(1:obj.maxConnections,:);
      
  end
  
end