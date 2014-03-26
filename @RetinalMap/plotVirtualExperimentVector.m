function fig = plotVirtualExperimentVector(obj,axisDir,titleStr)


  if(obj.plotFigures == 0)
    disp('plotVirtualInjectionExperimentVector: plotFigures = 0, hiding figures!')
    visFlag = 'off';
  else
    visFlag = 'on';
  end
  
  [SCdistance,flipFlag] = obj.virtualInjectionVectorExperiment(axisDir);
  
  % Now we want to bin by distance, to see what is probability to
  % flip as a function of distance between the two injections
  
  edges = linspace(0,1,25);
  Pflip = zeros(size(edges));
  
  for i = 1:numel(edges)-1
    idx = find(edges(i) <= SCdistance & SCdistance < edges(i+1));
    Pflip(i) = mean(flipFlag(idx));
  end
  
  % Plot the results
  fig = figure;
  stairs(edges,Pflip,'k-')
  switch(axisDir)
    case {'AP','ap'}
     xlabel('AP distance between injections')
    case {'ML','ml'}
     xlabel('ML distance between injections')
  end
     
  ylabel('Fraction of flips')
  
  
  fName = sprintf('%s/%s-virtual-exp-vector-%s.pdf', ...
                  obj.figurePath, obj.simName, axisDir);
  
  switch(obj.plotFigures)
    case {0,1}
      saveas(gcf,fName,'pdf');
    case 2
      obj.addParametersToFigure(strrep(fName,'.pdf','.eps'));
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
end
