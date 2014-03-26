function fig = plotGradients1D(obj)

  assert(strcmp(obj.eyeType,'disk'));

  if(obj.plotFigures == 0)
    disp('plotGradients1D: plotFigures = 0, hiding figures!')
    visFlag = 'off';
  else
    visFlag = 'on';
  end
  
  fig(1) = figure('visible',visFlag);
  
  if(~isempty(obj.RGCEphA))
    subplot(2,2,1)
    plot(obj.RGCnt, obj.RGCEphA, 'r.')
    xlabel('Temporal - Nasal')
    ylabel('EphA')
    title('Retina (forward)')
    set(gca,'xdir','reverse')
    box off
  else
    disp('No RGC EphA stored, ignoring')
  end

  if(~isempty(obj.RGCEphB))
    subplot(2,2,2)
    plot(obj.RGCdv, obj.RGCEphB, 'b.')
    xlabel('Ventral - Dorsal')
    ylabel('EphB')
    title('Retina (forward)')
    set(gca,'xdir','reverse')  
    box off
  else
    disp('No RGC EphB stored, ignoring')
  end
  
  if(obj.typeFlag > 1)
  
    if(~isempty(obj.RGCephrinA))
      subplot(2,2,3)
      plot(obj.RGCnt, obj.RGCephrinA, 'r.')
      xlabel('Temporal - Nasal')
      ylabel('ephrinA')
      title('Retina (reverse)')
      set(gca,'xdir','reverse')    
      box off
    else
      disp('No RGC ephrin-A stored, ignoring')
    end
    
    if(~isempty(obj.RGCephrinB))
      subplot(2,2,4)
      plot(obj.RGCdv, obj.RGCephrinB, 'b.')
      xlabel('Ventral - Dorsal')  
      ylabel('ephrinB')
      title('Retina (reverse)')
      set(gca,'xdir','reverse')    
      box off
    else
      disp('No RGC ephrin-B stored, ignoring')
    end
  end
  
  fig(2) = figure('visible',visFlag);
  if(~isempty(obj.SCephrinA))
    subplot(2,2,1)
    plot(obj.SCap, obj.SCephrinA, 'r.')
    xlabel('Anterior - Posterior')
    ylabel('ephrinA')
    title('SC (forward)')
    box off
  else
    disp('No SC ephrin-A stored, ignoring')
  end
  
  if(~isempty(obj.SCephrinB))
    subplot(2,2,2)
    plot(obj.SCml, obj.SCephrinB, 'b.')
    xlabel('Medial - Lateral')
    ylabel('ephrinB')
    title('SC (forward)')
    box off
  else
    disp('No SC ephrin-B stored, ignoring')
  end
  
  if(obj.typeFlag > 1)
  
    if(~isempty(obj.SCEphA))
      subplot(2,2,3)
      plot(obj.SCap, obj.SCEphA, 'r.')
      xlabel('Anterior - Posterior')
      ylabel('EphA')
      title('SC (reverse)')
      box off
    else
      disp('No SC EphA stored, ignoring')
    end
    
    if(~isempty(obj.SCEphB))
      subplot(2,2,4)
      plot(obj.SCml, obj.SCEphB, 'b.')
      xlabel('Medial - Lateral')  
      ylabel('EphB')  
      title('SC (reverse)')  
      box off
    else
      disp('No SC EphB stored, ignoring')
    end
  end
  
  fNameRet = sprintf('%s/%s-gradients-retina.pdf', ...
                     obj.figurePath, obj.simName);
  fNameSC = sprintf('%s/%s-gradients-SC.pdf', ...
                     obj.figurePath, obj.simName);

  switch(obj.plotFigures)
    case {0,1}
      saveas(fig(1),fNameRet,'pdf');
      saveas(fig(2),fNameSC,'pdf');      
    case 2
      figure(fig(1));
      obj.addParametersToFigure(strrep(fNameRet,'.pdf','.eps'));
      figure(fig(2));
      obj.addParametersToFigure(strrep(fNameSC,'.pdf','.eps'));      
  end
  
end