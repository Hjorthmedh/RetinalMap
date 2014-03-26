function fig = plotMapImage(obj,refImageName,interactiveMap,flipImg)

  if(~exist('refImageName') | isempty(refImageName))
    refImageName = 'images/grid-red-green-blue-yellow-black.png';
  end
  
  if(~exist('flipImg'))
    flipImg = 0;
  end

  if(~exist('interactiveMap'))
    interactiveMap = 0;
  end  
  
  [RGCcol,SCcol] = obj.getNeuronColoursFromImage(refImageName,flipImg);
  
  fig = obj.plotMapReverse(interactiveMap,SCcol,RGCcol);
  
  fName = sprintf('%s/%s-map-image-iter-%d.pdf', ...
                  obj.figurePath, obj.simName, obj.curStep/obj.nSC);

  switch(obj.plotFigures)
    case {0,1}
      saveas(gcf,fName,'pdf');
    case 2
      obj.addParametersToFigure(strrep(fName,'.pdf','.eps'));
  end
 
end