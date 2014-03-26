function loadExperimentConfigFile(obj,configFile)

  fprintf('Loading config file: %s\n', configFile)
  fid = fopen(configFile);
  
  str = fgets(fid);
  lineCtr = 0;
  
  while(str ~= -1)
    fprintf(str)
    try
      eval(str);
    catch e
      getReport(e)
      keyboard
    end
    str = fgets(fid);    
    lineCtr = lineCtr + 1;
  end
 
  fclose(fid);

  fprintf('Read %d lines\n', lineCtr)
  
end