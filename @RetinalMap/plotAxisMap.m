% This shows all synapses!
% For example: RGCX = RGCnt, SCY = SCap, nBins = 100
%
% SCrange allows us to limit the part of SC included in analysis
%
% axisRange = [1/3 2/3], uses middle third of retina

function plotAxisMap(obj, axisDir, nBins, axisRange, colorbarFlag, plotFlag)

  disp('plotAxisMap deprecated, please use plotAxisProjection')
  
  if(~exist('nBins') | isempty(nBins))
    nBins = 100; %200;
  end
  
  if(~exist('axisRange') | isempty(axisRange))
    axisRange = [0 1];
  end
 
  if(~exist('colorbarFlag') | isempty(colorbarFlag))
    colorbarFlag = 1;
  end
  
  if(~exist('plotFlag') | isempty(plotFlag))
    plotFlag = obj.plotFigures;
  end

  switch(axisDir)
    case {'NT','AP','nt','ap'}
      X = obj.RGCnt;
      Y = obj.SCap;
      xLabelName = 'Nasal-Temporal';
      yLabelName = 'Anterior-Posterior';
      xLabelAlt = {'N','T'};
      yLabelAlt = {'A','P'};
      Zfilt = obj.SCml; % If we only want to include part of SC
      Xmax = max(1,max(X));
      Ymax = max(1,max(Y));
      setAxis = true;
    case {'DV','ML','dv','ml'}
      X = obj.RGCdv;
      Y = obj.SCml;
      xLabelName = 'Dorsal-Ventral';
      yLabelName = 'Medial-Lateral';
      xLabelAlt = {'D','V'};
      yLabelAlt = {'M','L'};
      Zfilt = obj.SCap; % If we only want to include part of SC
      Xmax = max(1,max(X));
      Ymax = max(1,max(Y));
      setAxis = true;
    case {'EphA','ephrinA'}
      X = obj.RGCEphA;
      Y = obj.SCephrinA;
      xLabelName = 'Retinal EphA conc';
      yLabelName = 'SC ephrinA conc';
      xLabelAlt = {};
      yLabelAlt = {};
      Zfilt = obj.SCap;
      Xmax = max(1,max(X));
      Ymax = max(1,max(Y));
      setAxis = true;      
    case {'EphB','ephrinB'}
      X = obj.RGCEphB;
      Y = obj.SCephrinB;
      xLabelName = 'Retinal EphB conc';
      yLabelName = 'SC ephrinB conc';
      xLabelAlt = {};
      yLabelAlt = {};
      Zfilt = obj.SCml;
      Xmax = max(1,max(X));
      Ymax = max(1,max(Y));
      setAxis = true;      
    otherwise
      fprintf('Unknown axis %s, use AP,NT,DV or ML.\n', axisDir)
      return
  end

  if(Xmax == 0 | Ymax == 0)
    fprintf('plotAxisMap: Aborting since Xmax = %.2f, Ymax = %.2f\n', ...
            Xmax, Ymax);
    text(0.5,0.2,'One variable is constant')
    return
  end
  
  Zwidth = max(Zfilt(:)) - min(Zfilt(:));
  Zmin = min(Zfilt(:)) + axisRange(1)*Zwidth;
  Zmax = min(Zfilt(:)) + axisRange(2)*Zwidth;
  
  fprintf('Using SC axisRange from %d to %d (%d %d)\n', ...
          axisRange(1), axisRange(2), Zmin, Zmax)
  
  % colormap('hot');
  col = colormap('gray');
  colormap(col(end:-1:1,:));
  
  Xedges = linspace(0,Xmax,nBins);
  Yedges = linspace(0,Ymax,nBins);
  
  dx = Xedges(2)-Xedges(1);
  dy = Yedges(2)-Yedges(1);
  
  axisMap = zeros(nBins,nBins);

  filterCtr = 0;
  
  for i = 1:obj.nSC
    for j = 1:obj.numPresynapticConnections(i)
      idx = obj.presynapticConnections(j,i);
      w = obj.presynapticWeight(j,i);
   
      if(Zmin <= Zfilt(i) & Zfilt(i) <= Zmax)
        x = X(idx);
        xi = floor(x/dx) + 1;
        y = Y(i);
        yi = floor(y/dy) + 1;
        
        try
          axisMap(yi,xi) = axisMap(yi,xi) + w;
        catch e
          getReport(e)
          keyboard
        end
          
      else
        filterCtr = filterCtr + w;
      end
        
    end
  end

  fprintf('Filtered out %d synapses (weighted)\n', filterCtr)
  
  imagesc([0 Xmax],[0 Ymax],axisMap)
  xlabel(xLabelName,'fontsize',18)
  ylabel(yLabelName,'fontsize',18)
  set(gca,'ydir','normal','fontsize',14)

  if(setAxis)
    if(~isempty(xLabelAlt))
      a = axis;
      set(gca,'xtick',[0.1 0.9],'xticklabel',xLabelAlt);
      set(gca,'ytick',[0.1 0.9],'yticklabel',yLabelAlt);
      set(gca,'ticklength', [0 0]);     
      xlabel([])
      ylabel([])
      set(gca,'fontsize',30)
    else
    
      xLabelVal = 0:0.5:1;
      yLabelVal = 0:0.5:1;
      
      set(gca,'xtick',xLabelVal);
      set(gca,'ytick',yLabelVal);
      
    end
  end
    
  if(colorbarFlag)
    colorbar
  end
  
  axis equal
  axis tight
  box off
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  fName = sprintf('%s/%s-%s-axis.pdf', ...
                  obj.figurePath, obj.simName, ...
                  axisDir);

  if(exist('plotFlag') & plotFlag)
    switch(obj.plotFigures)
      case {0,1}
        saveas(gcf,fName,'pdf');
      case 2
        obj.addParametersToFigure(strrep(fName,'.pdf','.eps'));
    end
  end
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
end