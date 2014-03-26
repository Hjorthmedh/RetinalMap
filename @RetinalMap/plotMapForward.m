function fig = plotMapForward(obj,interactiveFlag,markIsl2)

  if(~exist('markIsl2'))
    markIsl2 = 0;
  end
  
  if(~exist('interactiveFlag') | isempty(interactiveFlag))
    interactiveFlag = 0;
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  assert(strcmp(obj.eyeType,'disk'))
  
  RGCcol = obj.getRGCcolour(); % Use function so we can handle
                               % surgery rotations etc.
  SCcol = zeros(obj.nSC,3);
  
  for i = 1:obj.nSC
    idx = obj.presynapticConnections(1:obj.numPresynapticConnections(i),i);
    w = obj.presynapticWeight(1:obj.numPresynapticConnections(i),i);
    
    for j = 1:numel(idx)
      try
        SCcol(i,:) = SCcol(i,:) + RGCcol(idx(j),:)*w(j);
      catch e
        getReport(e)
        keyboard
      end
    end
    
    SCcol(i,:) = SCcol(i,:) / sum(w);
    
  end
  
  fig = obj.plotMapReverse(interactiveFlag,SCcol,RGCcol,markIsl2);

  if(~exist(obj.figurePath))
    mkdir(obj.figurePath)
  end
  
  fName = sprintf('%s/%s-map-forward-iter-%d.pdf', ...
                  obj.figurePath, obj.simName, obj.curStep/obj.nSC);

  switch(obj.plotFigures)
    case {0,1}
      saveas(gcf,fName,'pdf');
    case 2
      obj.addParametersToFigure(strrep(fName,'.pdf','.eps'));
  end
    
end