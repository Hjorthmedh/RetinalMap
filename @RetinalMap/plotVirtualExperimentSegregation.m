% refCurves, one row per set of logistic parameters
% refColour, one row per set of colours
% lengthScale, one element per ref curve, for adult 2000 should be ok
%
% If legendText is given, then the first element in the cell array
% is the text for the model curve, the second to N+1 are for the N
% reference curves.

function [k,fig] = plotVirtualExperimentSegregation(obj,axisDir,titleText, ...
                                                    refCurves,refColour,lengthScale,legendText)

  nRep = 200;

  switch(axisDir)
    case {'NT','AP','nt','ap'}
      [d,sep,k,kLow,kHigh] = ...
          obj.virtualInjectionSegregationExperiment('MLfix',nRep);
    case {'DV','ML','dv','ml'}
      [d,sep,k,kLow,kHigh] = ...
          obj.virtualInjectionSegregationExperiment('APfix',nRep);
  end
  
  if(obj.plotFigures == 0)
    disp('plotVirtualExperimentSegregation: plotFigures = 0, hiding figures!')
    visFlag = 'off';
  else
    visFlag = 'on';
  end

  fig = figure('visible',visFlag);
  x = linspace(0,1,100);
  
  if(exist('refCurves') & ~isempty(refCurves))
    % Add reference curve(s) to the plot
    for i = 1:size(refCurves,1)
      if(exist('refColour') & ~isempty(refColour))
        colour = refColour(i,:);
      else
        colour = [1 1 1]*0.6;
      end
      
      if(exist('lengthScale') & ~isempty(lengthScale))
        scale = lengthScale(i);
      else
        scale = 1;
      end
      
      yRef = logistic(refCurves(i,:),x*scale);
      p(1+i) = plot(x,yRef,'-','color',colour,'linewidth',2);
      hold on
    end
    
  end
  
  plot(d,sep,'k.');
  hold on
  y = logistic(k,x);
  p(1) = plot(x,y,'k-','linewidth',2);
  hold off
  
  ylabel('Segregation in Retina','fontsize',24) 
  % Scale factor for 60 days old WT mouse, Dan Lyngholm personal
  % communication (email 18 july 2012). See fitSeggregationData.m
  % for entire table of sizes.
  title(sprintf('Segregation: inflection = %.3f (%.3f-%.3f), slope = %.3f (%.3f-%.3f)', ...
                k(1),kLow(1),kHigh(1),k(2),kLow(2),kHigh(2)))
  
  switch(axisDir)
    case {'NT','AP','nt','ap'}
      xlabel('Anterior - Posterior separation','fontsize',24)
    case {'DV','ML','dv','ml'}
      xlabel('Medial - Lateral separation','fontsize',24)
  end
  set(gca,'fontsize',20)
  
  box off
  
  if(exist('titleText') & ~isempty(titleText))
    title(titleText)
  end
  
  if(exist('legendText'))
    legend(p,legendText);
  end
  
  fName = sprintf('%s/%s-virtual-exp-segregation-%s.pdf', ...
                  obj.figurePath, obj.simName, axisDir);
  
  switch(obj.plotFigures)
    case {0,1}
      saveas(gcf,fName,'pdf');
    case 2
      obj.addParametersToFigure(strrep(fName,'.pdf','.eps'));
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function y = logistic(k,x)
    % y = -0.5./(1+(x/inflection).^slope)+1;
    y = 1 - 0.5./(1+(x/k(1)).^k(2));
        
  end  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
end
