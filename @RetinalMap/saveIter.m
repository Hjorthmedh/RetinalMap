function saveIter(obj,fileName)  

  if(~exist('fileName'))
    iterPath = sprintf('%s/ITER', obj.dataPath);

    fileName = sprintf('%s/%s-iter-%d.mat', ...
                       iterPath, obj.simName, obj.curStep);
    
    
    if(~exist(iterPath))
       fprintf('%s missing, creating.\n')
       mkdir(iterPath);
    end
  end
  
  obj.saveState(fileName);
  
end