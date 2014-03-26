function saveState(obj, fileName)
    
  if(~exist('fileName'))
    fileName = sprintf('%s/%s.mat', obj.dataPath, obj.simName);

    if(~exist(obj.dataPath))
      fprintf('%s missing, creating.\n', obj.dataPath)
      mkdir(obj.dataPath)
    end
  end

  old = struct('simName', []);
  
  allFields = fieldnames(obj); 
  
  % These can be recalculated easilly, no point storing
  excludeFields = {'UAct','CAct','neighbourRGC','neighbourSC'};
  
  saveFields = setdiff(allFields,excludeFields);
  
  for i = 1:numel(saveFields)
    eval(sprintf('old.%s = obj.%s;', saveFields{i}, saveFields{i}));
  end
  
  fprintf('Saving state to %s\n', fileName)

  if(obj.HDF5)
    save(fileName,'old','-v7.3'); % HDF5
  else
    save(fileName,'old');
  end
  
end