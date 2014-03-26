function run(obj)

  % Make sure the C and U tables exist!
  assert(~isempty(obj.CAct))
  assert(~isempty(obj.UAct))
  
  % Runs the full simulation
  tic
  if(obj.curStep == 0)
    obj.saveIter();
  end
  
  while(obj.curStep < obj.nSteps)
    iterSteps = min(obj.reportStep,obj.nSteps-obj.curStep);
    % obj.step(iterSteps);
    obj.stepMex(iterSteps);
    obj.saveIter();
  end
  toc

  obj.saveState();
  fprintf('Simulation %s done.\n', obj.simName)
  
end