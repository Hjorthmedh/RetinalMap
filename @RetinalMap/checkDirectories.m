function checkDirectories(obj)

  if(~exist(obj.dataPath))
    fprintf('Creating data directory: %s\n', obj.dataPath)
    mkdir(obj.dataPath)
  end
  
  if(~exist(obj.figurePath))
    fprintf('Creating figure directory: %s\n', obj.figurePath)
    mkdir(obj.figurePath)
  end

end