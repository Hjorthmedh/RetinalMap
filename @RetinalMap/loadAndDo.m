% Example usage:
% r = RetinalMap();
% r.loadAndDo('SAVE/ComparePhenotype','Comp*mat','plotAxisProjection(''NT'',false,false,''hist'',true); close all;')

function loadAndDo(obj,directory,pattern,action)

  files = dir(sprintf('%s/%s', directory,pattern));
  
  for i = 1:numel(files)
    fName = sprintf('%s/%s',directory,files(i).name);
    r = RetinalMap(fName);
    cmd = sprintf('r.%s;', action);
    eval(cmd);
  end

end