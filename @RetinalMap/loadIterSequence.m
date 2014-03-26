% This function returns a set of RetinalMap structs

function iterSequence = loadIterSequence(obj,skipTables,iter, ...
                                         dataPath,fileNameMask,uniqueFlag)

  simID = 0;
  
  if(~exist('fileNameMask'))
    fileNameMask = sprintf('%s-iter-*.mat',obj.simName);
    simID = obj.simID;
  end
  
  if(~exist('dataPath'))
    dataPath = obj.dataPath;
    simID = obj.simID;
  end
  
  if(~exist('skipTables') | isempty(skipTables))
    skipTables = 1;
  end
  
  if(~exist('uniqueFlag'))
    uniqueFlag = true;
  end
 
  if(strfind(fileNameMask,'%d') & ~isempty(iter))
    for i = 1:numel(iter)
      fName = sprintf(sprintf('%s/ITER/%s',dataPath,fileNameMask),iter(i));
      fprintf('Loading %s\n', fName)
      iterSequence(i) = RetinalMap();
      iterSequence(i).loadState(fName,skipTables);
    end
  else
    % Load files specified by the mask
    files = dir(sprintf('%s/ITER/%s',dataPath,fileNameMask));

    if(isempty(files))
      fprintf('No files found using mask: %s\n', fileNameMask)
      % keyboard
      iterSequence = [];
      return
    end
    
    for i = 1:numel(files)
      fName = sprintf('%s/ITER/%s',dataPath,files(i).name);
      fprintf('Loading %s\n', fName)
      iterSequence(i) = RetinalMap();
      iterSequence(i).loadState(fName,skipTables);
    end
  end
  
  % Sort the iterations
  curStep = cat(1,iterSequence.curStep);
  [~,idx] = sort(curStep);
  iterSequence = iterSequence(idx);

  if(simID == 0)
    simID = iterSequence(1).simID;
  end
    
  if(exist('iter') & ~isempty(iter))
    % Only return one specific iteration
    for i = 1:numel(iter)
      iterDist = abs(cat(1,iterSequence.curStep) - iter(i));
      [~,keepIdx(i)] = min(iterDist);
    end

    if(uniqueFlag)
      keepIdx = unique(keepIdx);
    end
    
    iterSequence = iterSequence(keepIdx);    
    fprintf('Only returning iteration %d\n', iterSequence.curStep)

    if(~skipTables)
      % We need to recalculate these tables
      iterSequence.updateNeighboursTables();
      iterSequence.calculateCU();
    end
  
  else
    allSimID = cat(1,iterSequence.simID);
    idx = find(allSimID - simID == 0);

    % Remove those that do not have the same simID as parent sim
    iterSequence = iterSequence(idx);

  end
  
  
end
