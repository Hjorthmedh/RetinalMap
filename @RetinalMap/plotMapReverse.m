% Normally SCcol and RGCcol are not specified by the user

function fig = plotMapReverse(obj,interactiveFlag,SCcol,RGCcol,markIsl2)

  silentFlag = 1;
  
  if(~exist('markIsl2'))
    markIsl2 = 0;
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % This allows the user to override the current colour selection
  if(~exist('SCcol'))
    
    % Calculate colour of each SC neuron
    SCb = 1-(obj.SCap - min(obj.SCap)) / (max(obj.SCap) - min(obj.SCap));
    SCg = 1-(obj.SCml - min(obj.SCml)) / (max(obj.SCml) - min(obj.SCml));

    SCcol = [zeros(size(SCb)) SCg SCb];
  end
  
  if(~exist('RGCcol'))
    RGCcol = zeros(obj.nRGC,3);

    for i = 1:obj.nSC
      for j = 1:obj.numPresynapticConnections(i)
        idx = obj.presynapticConnections(j,i);
        RGCcol(idx,:) = RGCcol(idx,:) + SCcol(i,:)*obj.presynapticWeight(j,i);
      end
    end
  
    for i=1:obj.nRGC
      RGCcol(i,:) = RGCcol(i,:) / double(obj.totalWeightRGC(i));
    end
  end
  
  if(obj.plotFigures == 0)
    disp('plotMapReverse: plotFigures = 0, hiding figures!')
    visFlag = 'off';
  else
    visFlag = 'on';
  end
  
  fig = figure('visible',visFlag);
  subplot(1,2,1);
  axesRetina = gca();
  subplot(1,2,2);
  axesSC = gca();
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if(exist('interactiveFlag') && interactiveFlag == 1 ...
     && strcmp(obj.eyeType,'disk'))
    
    set(fig,'windowbuttonmotionfcn',@mouseHandle);
    set(fig,'windowbuttondownfcn',@markPoint);
    
    % presynapticConnections show the connections SC -> RGC
    % Here we calculate RGC -> SC connectivity
    [mapR,mapWeightR,numConR] = RGCtoSCmap();
    
    % Centroid of the RGC neurons projecting to a particular SC
    [centroidRGCnt, centroidRGCdv] = calculateRGCcentroid();
    [maxRGCnt, maxRGCdv] = calculateRGCmax();
    
  else
    interactiveFlag = 0;
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  oldPointsRGC = [];
  oldPointsSC = [];
  oldRGCaxon = [];
  oldRGCcentroid = [];
  oldRGCmax = [];
  
  function mouseHandle(source, event)

    assert(strcmp(obj.eyeType,'disk'));
    
    if(~isempty(oldPointsRGC))
      try 
        delete(oldPointsRGC);
      catch,end

      oldPointsRGC = [];    
    end
    
    if(~isempty(oldPointsSC))
      try
        delete(oldPointsSC);
      catch,end
      oldPointsSC = [];
    end

    if(~isempty(oldRGCcentroid))
      try
        delete(oldRGCcentroid);
      catch, end
      oldRGCcentroid = [];
    end
    
    if(~isempty(oldRGCmax))
      try
        delete(oldRGCmax);
      catch, end
      oldRGCmax = [];
    end    
    
    if(~isempty(oldRGCaxon))
      try
      delete(oldRGCaxon);
      catch, end
      oldRGCaxon = [];
    end
    
    figure(fig)
    
    set(fig,'currentaxes',axesRetina);
    xy = get(axesRetina,'currentpoint');
    x = xy(1,1);
    y = xy(1,2);
    
    a = axis();
    
    if(a(1) <= x & x <= a(2) ...
       & a(3) <= y & y <= a(4))
      % Valid point in the retina, mark it

      [idx,dist] = findClosestRGC(x,y);

      if(dist > 0.1)
        idx = [];
      else
        if(~silentFlag)
          fprintf('Showing RGC %d and its postsynaptic partners\n', idx)
        end
      end     
      

      
      % plotRetinalMap();
      
      set(fig,'currentaxes',axesRetina);
      oldPointsRGC = plot(obj.RGCnt(idx),obj.RGCdv(idx), ...
                          'wo','markersize',8);     
     
      % plotSCmap();

      n = numConR(idx);
      
      oldPointsSC = zeros(n,1);
      set(fig,'currentaxes',axesSC);

      for i = 1:n
        jdx = mapR(i,idx);
        
        oldPointsSC(i) = plot(obj.SCap(jdx), obj.SCml(jdx), ...
                              'wo', 'markersize', mapWeightR(i,idx)+4);
      end
      
      % Also mark where the axon enters
      
      if(~isempty(idx))
        if(~isempty(obj.RGCml))
          oldRGCaxon = plot([min(obj.SCap(:)) max(obj.SCap(:))], ...
                            obj.RGCml(idx)*[1 1], 'k-');
        end
      end
    end

   
    set(fig,'currentaxes',axesSC);
    xy = get(axesSC,'currentpoint');
    x = xy(1,1);
    y = xy(1,2);  
 
    a = axis();

    if(a(1) <= x & x <= a(2) ...
       & a(3) <= y & y <= a(4))
      % Valid point in the SC, mark it
  
      [idx,dist] = findClosestSC(x,y);
      
      if(dist > 0.1)
        idx = [];
      else
        if(~silentFlag)
          fprintf('Showing SC %d and its presynaptic partners\n', idx)
        end
      end
      

      
      % plotSCmap();

      set(fig,'currentaxes',axesSC);

      oldPointsSC = plot(obj.SCap(idx),obj.SCml(idx),'wo');
      
      % plotRetinalMap();

      n = obj.numPresynapticConnections(idx);      
      oldPointsRGC = zeros(n,1);

      set(fig,'currentaxes',axesRetina);      
      for i = 1:n
        jdx = obj.presynapticConnections(i,idx);
        oldPointsRGC(i) = plot(obj.RGCnt(jdx), obj.RGCdv(jdx), ...
                               'wo', 'markersize', obj.presynapticWeight(i,idx)+4);
        
      end
      % plotCentroid
      oldRGCcentroid = plot(centroidRGCnt(idx), centroidRGCdv(idx), ...
                            'w*', 'markersize',10);

      oldRGCmax = plot(maxRGCnt(idx), maxRGCdv(idx), ...
                       'y*', 'markersize',10);
      
    end
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  nMarked = 0;
  markColor = ['rg'];
  
  function markPoint(source, event)
    
    markCol = markColor(mod(nMarked,numel(markColor))+1);
    
    % This allows the user to mark a neuron
    
    set(fig,'currentaxes',axesRetina);
    xy = get(axesRetina,'currentpoint');
    x = xy(1,1);
    y = xy(1,2);
    
    a = axis();
    
    if(a(1) <= x & x <= a(2) ...
       & a(3) <= y & y <= a(4))
      % Valid point in the retina, mark it

      [idx,dist] = findClosestRGC(x,y);

      if(dist > 0.1)
        idx = [];
      else
        fprintf('Marking RGC %d and its postsynaptic partners\n', idx)
        nMarked = nMarked + 1;
      end     

      set(fig,'currentaxes',axesRetina);
      plot(obj.RGCnt(idx),obj.RGCdv(idx), [markCol 'o'],'markersize',8);     

      set(fig,'currentaxes',axesSC);
      n = numConR(idx);

      for i = 1:n
        jdx = mapR(i,idx);
        
        plot(obj.SCap(jdx), obj.SCml(jdx), ...
             [markCol 'o'], 'markersize', mapWeightR(i,idx)+4);
      end
  
      % Also mark where the axon enters
      
      if(~isempty(idx) & ~isempty(obj.RGCml))
        plot([min(obj.SCap(:)) max(obj.SCap(:))], ...
             obj.RGCml(idx)*[1 1], [markCol '-']);
      end
      
    end

    
    set(fig,'currentaxes',axesSC);
    xy = get(axesSC,'currentpoint');
    x = xy(1,1);
    y = xy(1,2);  
 
    a = axis();

    if(a(1) <= x & x <= a(2) ...
       & a(3) <= y & y <= a(4))
      % Valid point in the SC, mark it
  
      [idx,dist] = findClosestSC(x,y);
      
      if(dist > 0.1)
        idx = [];
      else
        fprintf('Marking SC %d and its presynaptic partners\n', idx)
        nMarked = nMarked + 1;
      end    
    
      set(fig,'currentaxes',axesSC);
      plot(obj.SCap(idx),obj.SCml(idx), [markCol 'o']);

      n = obj.numPresynapticConnections(idx);      
      
      set(fig,'currentaxes',axesRetina);      
      for i = 1:n
        jdx = obj.presynapticConnections(i,idx);
        plot(obj.RGCnt(jdx), obj.RGCdv(jdx), ...
             [markCol 'o'], 'markersize', obj.presynapticWeight(i,idx)+4);
        
      end
      % plotCentroid
      plot(centroidRGCnt(idx), centroidRGCdv(idx), ...
           [markCol '*'], 'markersize',10);
      
    end
     
      
      
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function [idx,dist] = findClosestRGC(x,y)
    
    switch(obj.eyeType)
      
      case 'sphere'
        disp('Sphere interactive plots not currently implemented')
        keyboard
        
      case 'disk'
        [dist,idx] = min(sqrt((obj.RGCnt - x).^2 + (obj.RGCdv - y).^2));
        
      otherwise
        fprintf('Unknown eye type: %s\n', obj.eyeType)
        keyboard
        
    end
          
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function [idx,dist] = findClosestSC(x,y)
    [dist,idx] = min(sqrt((obj.SCap - x).^2 + (obj.SCml - y).^2));
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function [map,mapWeight,numCon] = RGCtoSCmap()
    
    map = zeros(ceil(max(obj.totalWeightRGC)),obj.nRGC);
    mapWeight = zeros(ceil(max(obj.totalWeightRGC)),obj.nRGC);
    targetCtr = zeros(1,obj.nRGC);
    
    for i = 1:obj.nSC
      for j = 1:obj.numPresynapticConnections(i)
        rIdx = obj.presynapticConnections(j,i);
        rW = obj.presynapticWeight(j,i);
        
        targetCtr(rIdx) = targetCtr(rIdx) + 1;
        
        map(targetCtr(rIdx),rIdx) = i;
        mapWeight(targetCtr(rIdx),rIdx) = rW;
        
      end
      
    end
    
    numCon = targetCtr;
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function plotRetinalMap()
  
    set(fig,'currentaxes',axesRetina);
    
    switch(obj.eyeType)
      case 'sphere'
        hold off
        polar(linspace(0,2*pi,100),obj.maxTheta*ones(1,100),'r-')
        hold on
  
        for i = 1:obj.nRGC
          r = obj.RGCtheta(i);
          v = obj.RGCphi(i);
        
          if(any(isnan(RGCcol(i,:))))
            plot(r*cos(v),r*sin(v), 'ko');
          else
            plot(r*cos(v),r*sin(v), '.', 'color', RGCcol(i,:), ...
                 'markersize',12);
          end
        end

        title('Remember to add axes labels')
      
        axis equal         
      
      case 'disk'
        hold off
        if(markIsl2 & ~isempty(obj.Isl2PositiveRGC))
          plot(obj.RGCnt(obj.Isl2PositiveRGC), ...
               obj.RGCdv(obj.Isl2PositiveRGC), ...
               'o','color',[1 0 0]*0.8, ...
               'markersize', 6);
          hold on
        end
        for i = 1:obj.nRGC
          if(any(isnan(RGCcol(i,:))))
            % plot(obj.RGCnt(i), obj.RGCdv(i), 'ko');
            plot(obj.RGCnt(i), obj.RGCdv(i), ...
                 '.','markersize',6,'color',[1 1 1]*0.6);
          else
            plot(obj.RGCnt(i), obj.RGCdv(i),...
                 '.', 'color', RGCcol(i,:), ...
                 'markersize', 12);
          end
          hold on        
        end

        % We want T-N on X axis, and V-D on Y-axis, need to reverse them
        set(gca, 'xdir','reverse','ydir','reverse');        
        
        %ylabel('Ventral - Dorsal','fontsize',16)
        %xlabel('Temporal - Nasal','fontsize',16)
        set(gca,'fontsize',30)
    
        axis equal    
        mx = max(max(obj.RGCnt),1);
        my = max(max(obj.RGCdv),1);
        axis([0 mx 0 my])
        
        set(gca,'xtick',[0.1 0.9]*mx,'xticklabel',{'N','T'});
        set(gca,'ytick',[0.1 0.9]*my,'yticklabel',{'D','V'});
        set(gca,'ticklength', [0 0]);        
        box off
            
      otherwise
        fprintf('plotMapReverse: Unknown eye type %s\n', obj.eyeType)
        keyboard  
    end
  
    if(interactiveFlag)
      set(gca,'color',0.2*[1 1 1]);
    end    
    
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function plotSCmap()

    set(fig,'currentaxes',axesSC);    

    hold off
    for i = 1:numel(obj.SCml)
      if(obj.numPresynapticConnections(i) > 0)
        plot(obj.SCap(i),obj.SCml(i),'.', ...
             'color',SCcol(i,:), ...
             'markersize',12);
      else
        % plot(obj.SCap(i),obj.SCml(i),'ko');     
        plot(obj.SCap(i),obj.SCml(i), ...
             '.','markersize',6,'color',[1 1 1]*0.6);
      end
      hold on
    end
  

    %xlabel('Anterior - Posterior','fontsize',16)
    %ylabel('Medial - Lateral','fontsize',16)
    set(gca,'fontsize',30)
  
    if(0)
      % Mark the axons
      
      for i = 1:numel(obj.RGCml)
        plot([0 1],obj.RGCml(i)*[1 1],'k-')
      end
    end
    
    axis equal
    mx = max(max(obj.SCap),1);
    my = max(max(obj.SCml),1);
    axis([0 mx 0 my])
    set(gca,'xtick',[0.1 0.9]*mx,'xticklabel',{'A','P'});
    set(gca,'ytick',[0.1 0.9]*my,'yticklabel',{'M','L'});
    set(gca,'ticklength', [0 0]);
    box off
    
    if(interactiveFlag)
      set(gca,'color',0.2*[1 1 1]);
    end
  
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function [centroidRGCnt, centroidRGCdv] = calculateRGCcentroid()
    centroidRGCnt = zeros(obj.nSC,1);
    centroidRGCdv = zeros(obj.nSC,1);

    for i = 1:obj.nSC
      totalWeightSC = ...
          sum(obj.presynapticWeight(1:obj.numPresynapticConnections(i),i));

      for j = 1:obj.numPresynapticConnections(i)
        centroidRGCnt(i) = centroidRGCnt(i) ...
            + obj.presynapticWeight(j,i) ...
              * obj.RGCnt(obj.presynapticConnections(j,i));
        centroidRGCdv(i) = centroidRGCdv(i) ...
            + obj.presynapticWeight(j,i) ...
              * obj.RGCdv(obj.presynapticConnections(j,i));
      end
    
      centroidRGCnt(i) = centroidRGCnt(i) / totalWeightSC;
      centroidRGCdv(i) = centroidRGCdv(i) / totalWeightSC;
      
    end
    
  end

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function [maxRGCnt, maxRGCdv] = calculateRGCmax()
    maxRGCnt = zeros(obj.nSC,1);
    maxRGCdv = zeros(obj.nSC,1);

    for i = 1:obj.nSC
      totalWeightSC = ...
          sum(obj.presynapticWeight(1:obj.numPresynapticConnections(i),i));

      nPre = obj.numPresynapticConnections(i);
      [~, idx] = max(obj.presynapticWeight(1:nPre,i));
      if(~isempty(idx))
        maxRGCnt(i) = obj.RGCnt(obj.presynapticConnections(idx,i));
        maxRGCdv(i) = obj.RGCdv(obj.presynapticConnections(idx,i));
      else
        maxRGCnt(i) = NaN;
        maxRGCdv(i) = NaN;
      end
    end
    
  end
      
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  plotRetinalMap();
  %title(sprintf('Iteration %d', obj.curStep/obj.nSC))
  plotSCmap();
  
end