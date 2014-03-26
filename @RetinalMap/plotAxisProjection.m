function fig = plotAxisProjection(obj,showAxis,useAll,reuseFigFlag,plotType,saveFigFlag,weightScale,markerAlpha,showExpFlag)

  if(~exist('showAxis') | isempty(showAxis))
    showAxis = 'NT';
  end
  
  if(~exist('useAll') | isempty(useAll))
    useAll = true;
  end
  
  if(~exist('reuseFigFlag') | isempty(reuseFigFlag))
    reuseFigFlag = false;
  end
    
  if(~exist('weightScale') | isempty(weightScale))
    weightScale = 1.0;
  end
  
  if(~exist('markerAlpha') | isempty(markerAlpha))
    markerAlpha = 0.05;
  end
  
  if(~exist('saveFigFlag') | isempty(saveFigFlag))
    saveFigFlag = false;
  end
  
  if(~reuseFigFlag)
    if(obj.plotFigures)
      fig = figure;
    else
      fig = figure('visible','off');
    end
  else
    fig = gcf;
  end
    
  if(~exist('plotType') | isempty(plotType))
    plotType = 'alpha';
  end
  
  if(~exist('showExpFlag') | isempty(showExpFlag))
    showExpFlag = true;
  end

  if(useAll)
    disp('Using all RGC neurons')
  else
    disp('Restricting to middle third of the retina')
  end
  
  normalColour = [0 0 0];
  altColour = [0 0 1];
  
  switch(showAxis)
    case {'NT','nt','AP','ap'}
      X = obj.RGCnt;
      Y = obj.SCap;
      xLab = {'N', 'T'}; % 'Nasal - Temporal';
      yLab = {'A', 'P'}; %'Anterior - Posterior';
      axDir = 'NT';
      
      if(useAll)
        useIdx = 1:obj.nRGC;
      else
        useIdx = find(1/3 < obj.RGCdv & obj.RGCdv < 2/3);
      end
      
    case {'DV','dv','ML','ml'}
      X = obj.RGCdv;
      Y = obj.SCml;
      xLab = {'D','V'}; %'Dorsal - Ventral';
      yLab = {'M','L'}; %'Medial - Lateral';
      axDir = 'DV';
      
      if(useAll)
        useIdx = 1:obj.nRGC;
      else
        useIdx = find(1/3 < obj.RGCnt & obj.RGCnt < 2/3);
      end
      
    otherwise
      fprintf('Unknown axis: %s, use NT or DV\n')
      return
  end

  obj.convertConnectionTables('pre2post');

  nChar = numel(obj.phenotype);
  str = obj.phenotype(1:min(4,nChar));
  
  switch(str)
    case {'Isl2'}
      
      switch(plotType)
        case 'alpha'
          dotPlotAlpha(obj.Isl2PositiveRGC,useIdx,markerAlpha);
      
        case 'dot'
          colourP = altColour;
          colourM = normalColour;
        
          pIdx = intersect(obj.Isl2PositiveRGC,useIdx);
          mIdx = intersect(setdiff(1:obj.nRGC,obj.Isl2PositiveRGC),useIdx);
          dotPlotNoAlpha(pIdx,colourP);
          dotPlotNoAlpha(mIdx,colourM);
        
        case 'hist'
          histPlot(obj.Isl2PositiveRGC,useIdx);
          
        otherwise
          fprintf('Unknown plot type: %s (use alpha, dot or hist)\n', plotType)
          keyboard
      end
      
    otherwise
      
      switch(plotType)
        
        case 'alpha'
          % All normal cells, no alt
          dotPlotAlpha([],useIdx,markerAlpha);
        
        case 'dot'
          colour = normalColour;
          dotPlotNoAlpha(useIdx,colour,markerAlpha)
          
        case 'hist'
          histPlot([],useIdx);
          
        otherwise
          fprintf('Unknown plot type: %s (use alpha, dot or hist)\n', plotType)
          keyboard
          
      end
  end
  
  if(showExpFlag & strcmp(axDir,'NT'))
    switch(obj.phenotype)
      case 'WT'
        expData = load('LemkeIsl2Injections/wt.csv');
        expData = expData / 100; % Transform to our coordinate system
        %plot(expData(:,1),expData(:,2), ...
        %     '.','markersize',30,'color',[51 135 38]/255)
        plotExpData(expData(:,1),expData(:,2),1,[1 0 0])
        
      case {'Isl2homozygous','Isl2hom'}
        expData = load('LemkeIsl2Injections/ki.ki.csv');
        expData = expData / 100;
        
        plotExpData(expData(:,1),expData(:,3),1,[1 0 0]) % Isl2-
        plotExpData(expData(:,1),expData(:,2),1,[1 0 0]) % Isl2+
        
      case {'Isl2heterozygous','Isl2het'}
        expData = load('LemkeIsl2Injections/ki.+.csv');
        expData = expData / 100;
        
        plotExpData(expData(:,1),expData(:,3),1,[1 0 0]) % Isl2-
        plotExpData(expData(:,1),expData(:,2),1,[1 0 0]) % Isl2+
        plotExpData(expData(:,1),expData(:,4),1,[1 0 0]) % Unknown

      case {'Isl2homEphA4mm'}        
        
        expData = load('LemkeIsl2Injections/ki.ki_A4-.-.csv');
        expData = expData / 100;
        
        plotExpData(expData(:,1),expData(:,3),1,[1 0 0])
        plotExpData(expData(:,1),expData(:,2),1,[1 0 0])
                
      case {'Isl2homEphA4pm'}        
        
        expData = load('LemkeIsl2Injections/ki.ki_A4+.-.csv');
        expData = expData / 100;
        
        plotExpData(expData(:,1),expData(:,3),1,[1 0 0])
        plotExpData(expData(:,1),expData(:,2),1,[1 0 0])        
        
      case {'Isl2hetEphA4mm'}
        
        expData = load('LemkeIsl2Injections/ki.+_A4-.-.csv');
        expData = expData / 100;
        
        plotExpData(expData(:,1),expData(:,3),1,[1 0 0])
        plotExpData(expData(:,1),expData(:,2),1,[1 0 0])        
        
      case {'Isl2hetEphA4pm'}        
        
        expData = load('LemkeIsl2Injections/ki.+_A4+.-.csv');
        expData = expData / 100;
        
        plotExpData(expData(:,1),expData(:,3),1,[1 0 0])
        plotExpData(expData(:,1),expData(:,2),1,[1 0 0])
      
      case {'Isl2homEphA5mm'}        
        
      case {'Isl2homEphA5pm'}        
        
      case {'Isl2hetEphA5mm'}
        
      case {'Isl2hetEphA5pm'}        

        
        
    end
  end
  
  set(gca,'ydir','normal','fontsize',30)
     
  set(gca,'xtick',[0.1 0.9],'xticklabel',xLab)
  set(gca,'ytick',[0.1 0.9],'yticklabel',yLab)
  set(gca,'ticklength', [0 0]);
  
  if(useAll)
    title('All RGC included')
  end
  
  %xlabel(xLab,'fontsize',30);
  %ylabel(yLab,'fontsize',30);

  axis equal
  axis([0 1 0 1])
  
  box off
  
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if(~reuseFigFlag)
    if(saveFigFlag)
      % Plot the figure
  
      
      switch(obj.plotFigures)
        case {0,1}
          switch(plotType)
            case {'alpha'}
              fName = sprintf('%s/%s-%s-projection.png', ...
                              obj.figurePath, obj.simName, showAxis);
      
              fprintf('Saving to %s\n', fName)
              
              print(gcf,'-dpng',fName,'-r600');
            otherwise
              fName = sprintf('%s/%s-%s-projection.pdf', ...
                              obj.figurePath, obj.simName, showAxis);
      
              fprintf('Saving to %s\n', fName)

              % saveas(gcf,fName,'pdf')
              print(gcf,'-dpdf',fName,'-painters','-r1200');          
          end

        case 2
          obj.addParametersToFigure(strrep(fName,'.png','.eps'));
      end
    end
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function dotPlotNoAlpha(RGCidx,colour)

    hold on

    xAll = [];
    yAll = [];
    wAll = [];
    
    for i = 1:numel(RGCidx)
      rIdx = RGCidx(i);
      
      n = obj.numPostsynapticConnections(rIdx);
      idx = obj.postsynapticConnections(1:n,rIdx);
      w = obj.postsynapticWeight(1:n,rIdx);

      x = X(rIdx)*ones(n,1);
      y = Y(idx);
      
      xAll = [xAll; x];
      yAll = [yAll; y];
      wAll = [wAll; w];
    end
    
    scatter(xAll,yAll,10*wAll*weightScale,colour)
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function dotPlotAlpha(RGCidxAlt,includeIdx,markerAlpha)

    hold on
    r = 0.005*weightScale;
    v = linspace(0,2*pi,10);
    dx = r*cos(v);
    dy = r*sin(v);

    % Only loop over the RGC that are included
    for iFilt = 1:numel(includeIdx)
      
      rIdx = includeIdx(iFilt);
      
      n = obj.numPostsynapticConnections(rIdx);
      idx = obj.postsynapticConnections(1:n,rIdx);
      w = obj.postsynapticWeight(1:n,rIdx);

      x = X(rIdx)*ones(n,1);
      y = Y(idx);
      
      ws = sqrt(w); % We want area of marker to be proportional to weight
        
      if(ismember(rIdx,RGCidxAlt))
        % Special cells (if they exist) marked in red, e.g. Isl2+
        colour = altColour;
      else
        % Normal cells marked in black
        colour = normalColour;
      end
      
      for j = 1:n
        p = patch(x(j)+ws(j)*dx,y(j)+ws(j)*dy, colour, ...
                  'facecolor', colour, ...
                  'edgecolor', colour, ...
                  'alphadatamapping','direct', ...
                  'facealpha',markerAlpha, ...
                  'edgealpha', markerAlpha);
        set(p,'zdata',zeros(size(get(p,'xdata'))));
      end
    end
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function histPlot(RGCidxAlt, includeIdx)
  
    blackIdx = setdiff(includeIdx,RGCidxAlt);
    blueIdx = intersect(includeIdx, RGCidxAlt);

    hBlack = histPlotHelper(blackIdx);
    hBlue = histPlotHelper(blueIdx);
    
    % Default colour is white
    H = zeros(size(hBlack,1),size(hBlack,2),3);
    H(:,:,1) = - hBlack - hBlue;
    H(:,:,2) = - hBlack - hBlue;
    H(:,:,3) = - hBlack;    
    
    % Scale colours to range 0 to 1
    hMax = max(H(:));
    hMin = min(H(:));
    hRange = hMax - hMin;
    
    H = (H - hMin) / hRange;
    
    image([min(X) max(X)],[min(Y) max(Y)], H);
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function h = histPlotHelper(RGCidx)
    
    nBins = 100;
    
    h = zeros(nBins,nBins);
    ex = linspace(0,max(X),nBins);
    ey = linspace(0,max(Y),nBins);
    
    dx = ex(2)-ex(1);
    dy = ey(2)-ey(1);
    
    for iR = 1:numel(RGCidx)
      rIdx = RGCidx(iR);
      
      % Max so that if X = 0, we get index 1.      
      hx = max(1,ceil(X(rIdx)/dx));
      
      for in = 1:obj.numPostsynapticConnections(rIdx)
        idx = obj.postsynapticConnections(in,rIdx);
        w = obj.postsynapticWeight(in,rIdx);
        hy = max(1,ceil(Y(idx)/dy));
        
        try
          h(hy,hx) = h(hy,hx) + w;
        catch e
          getReport(e)
          keyboard
        end
      end
      
    end
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  
  function plotExpData(NT,AP,expAlpha,expColour,v)

    if(~exist('expAlpha'))
      expAlpha = 1;
    end
    
    if(~exist('expColour'))
      expColour = [51 135 38]/255;
    end
    
    
    
    r = 0.015;
    if(~exist('v') | isempty(v))
      v = linspace(0,2*pi,50);
    end
    dx = r*cos(v);
    dy = r*sin(v);
    
    for i = 1:numel(NT)

      p = patch(NT(i) + dx, AP(i) + dy, expColour, ...
                'facecolor', expColour, ...
                'edgecolor', [1 1 1], ...
                'alphadatamapping','none', ...
                'facealpha', expAlpha, ...
                'edgealpha', expAlpha);
      
      set(p,'zdata',i+0.1*ones(size(get(p,'xdata'))));      
    end
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
end