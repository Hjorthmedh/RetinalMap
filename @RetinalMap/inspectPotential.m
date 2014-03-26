% Plots all the RGC and SC neurons, and allow the user to click on
% one of them to see the potential.

function inspectPotential(obj,addFlag, SCidx)

  if(isempty(obj.UAct))
    disp('You need to have C and U tables')
    return
  end
  
    % Whether we do synapse addition or removal energy
  if(~exist('addFlag'))
    addFlag = 1;
  end

  if(exist('SCidx'))
    selectIdx = SCidx;
  else
    selectIdx = [];
  end
  
  disp('Currently only works with eye-disks')
  assert(strcmp(obj.eyeType,'disk'));

  fig = figure();
  %set(fig,'windowbuttonmotionfcn',@mouseHandle,'interruptible','off');
  set(fig,'windowbuttondownfcn',@mouseHandle,'interruptible','off');
  subplot(2,3,1)
  axesRetTot = gca(); 
  subplot(2,3,2);
  axesRetChem = gca();
  subplot(2,3,3);
  axesSC = gca();
  subplot(2,3,4);
  axesRetAct = gca();
  subplot(2,3,5);
  axesRetComp = gca();
  subplot(2,3,6);
  axesRetMap = gca();

  if(obj.nSC >= 5000)
    markSize = 12;
  else
    markSize = 18;
  end
  
  [RGCcol,SCcol] = getMapColor();
  
  plotFigure();

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function [RGCcol, SCcol] = getMapColor()
    RGCb = (obj.RGCnt - min(obj.RGCnt)) / (max(obj.RGCnt) - min(obj.RGCnt));
    RGCg = (obj.RGCdv - min(obj.RGCdv)) / (max(obj.RGCdv) - min(obj.RGCdv));
  
    RGCcol = [zeros(size(RGCb)) RGCg RGCb];
    SCcol = zeros(obj.nSC,3);
  
    for i = 1:obj.nSC
      
      if(obj.numPresynapticConnections(i) > 0)
        idx = obj.presynapticConnections(1:obj.numPresynapticConnections(i),i);
        w = obj.presynapticWeight(1:obj.numPresynapticConnections(i),i);
        
        for j = 1:numel(idx)
          SCcol(i,:) = SCcol(i,:) + RGCcol(idx(j),:)*w(j);
        end
    
        SCcol(i,:) = SCcol(i,:) / sum(w);

      else
        SCcol(i,:) = 0.7;
      end
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function plotFigure()

    figure(fig);
    
    if(~isempty(selectIdx))
      fprintf('Calculating for SC %d at ml = %.3f, ap = %.3f...', ... 
              selectIdx, obj.SCml(selectIdx), obj.SCap(selectIdx))
    end
      
    set(fig,'currentaxes',axesSC);
    cla, hold on
    % Plot all SC
    for i = 1:obj.nSC
      plot(obj.SCap(i),obj.SCml(i), ...
           '.','color',SCcol(i,:), 'markersize',markSize);
    end
       
    if(~isempty(selectIdx))
      plot(obj.SCap(selectIdx),obj.SCml(selectIdx), ...
           'wo', 'markersize', markSize);
    
      plot(obj.SCap(selectIdx),obj.SCml(selectIdx), ...
           'ro', 'markersize', markSize*0.8);
      
      % Lets just mark the neuron in big colours
      plot(0,0,'.','color',SCcol(selectIdx,:),'markersize', 40);
    end
    
    if(~isempty(selectIdx))
      if(addFlag)
        title('Adding a synapse')
      else
        title('Removing a synapse')
      end
    else
      title('SC map')
    end

    xlabel('Anterior-Posterior')    
    ylabel('Medial-Lateral')

    axis equal
    box off
    hold off
    
    if(~isempty(selectIdx))
    
      Etot = zeros(obj.nRGC,1);
      Echem = zeros(obj.nRGC,1);
      Eact = zeros(obj.nRGC,1);
      Ecomp = zeros(obj.nRGC,1);

      totCol = zeros(obj.nRGC,3);
      chemCol = zeros(obj.nRGC,3);
      actCol = zeros(obj.nRGC,3);
      compCol = zeros(obj.nRGC,3);   
    
      if(addFlag)
        % Energy to add a synapse
        for iRGC=1:obj.nRGC
          if(mod(iRGC,100) == 0)
            fprintf(' %d', iRGC)
          end
          try
            [Etot(iRGC),Echem(iRGC),Eact(iRGC),Ecomp(iRGC)] = ...
                obj.calculateSynapseAdditionEnergy(selectIdx,iRGC);
          catch e
            disp('In inspectPotential.m')
            getReport(e)
            keyboard
          end
        end      
      else
        % Energy to remove a synapse
        % We only need to loop through the RGC that are connected
        
        Etot = NaN*Etot;
        Echem = NaN*Echem;
        Eact = NaN*Eact;
        Ecomp = NaN*Ecomp;
        
        for iN=1:obj.numPresynapticConnections(selectIdx)
          iRGC = obj.presynapticConnections(iN,selectIdx);
          [Etot(iRGC),Echem(iRGC),Eact(iRGC),Ecomp(iRGC)] = ...
              obj.calculateSynapseRemovalEnergy(selectIdx,iRGC);
          
        end
               
      end
      
      % Find a scaling factor for the plots
      eMax = max([max(Etot),max(Echem),max(Eact),max(Ecomp)]);
      eMin = min([min(Etot),min(Echem),min(Eact),min(Ecomp)]);
      
      colormap('jet');
      colMap = colormap('jet');
      
      eSpace = linspace(eMin,eMax,64);

      if(addFlag)
        for i = 1:3
          totCol(:,i) = interp1(eSpace,colMap(:,i),Etot);
          chemCol(:,i) = interp1(eSpace,colMap(:,i),Echem);
          actCol(:,i) = interp1(eSpace,colMap(:,i),Eact);
          compCol(:,i) = interp1(eSpace,colMap(:,i),Ecomp);     
        end
      else
        totCol(:,:) = 0.8;
        chemCol(:,:) = 0.8;
        actCol(:,:) = 0.8;
        compCol(:,:) = 0.8;        

        idx = find(~isnan(Etot));

        for i = 1:3
          totCol(idx,i) = interp1(eSpace,colMap(:,i),Etot(idx));
          chemCol(idx,i) = interp1(eSpace,colMap(:,i),Echem(idx));
          actCol(idx,i) = interp1(eSpace,colMap(:,i),Eact(idx));
          compCol(idx,i) = interp1(eSpace,colMap(:,i),Ecomp(idx));     
        end
        
      end

      if(addFlag)
      
        set(fig,'currentaxes', axesRetTot);      
        cla, hold on
        plotEnergy(obj.RGCnt,obj.RGCdv,totCol,markSize);
        hold off
        title('Total energy change')
        ylabel('Ventral-Dorsal')
        xlabel('Temporal-Nasal')
        set(gca,'xdir','reverse','ydir','reverse')
        axis equal
        box off
        
        set(fig,'currentaxes', axesRetChem);
        cla, hold on
        plotEnergy(obj.RGCnt,obj.RGCdv,chemCol,markSize);
        hold off
        title('Chemical energy change')
        ylabel('Ventral-Dorsal')
        xlabel('Temporal-Nasal')
        set(gca,'xdir','reverse','ydir','reverse')
        axis equal
        box off
        
        set(fig,'currentaxes', axesRetAct);    
        cla, hold on
        plotEnergy(obj.RGCnt,obj.RGCdv,actCol,markSize);
        hold off
        title('Activity energy change')
        ylabel('Ventral-Dorsal')
        xlabel('Temporal-Nasal')
        set(gca,'xdir','reverse','ydir','reverse')
        axis equal
        box off
      
        set(fig,'currentaxes', axesRetComp);
        cla, hold on
        plotEnergy(obj.RGCnt,obj.RGCdv,compCol,markSize);
        hold off
        title('Competition energy change')
        ylabel('Ventral-Dorsal')
        xlabel('Temporal-Nasal')
        set(gca,'xdir','reverse','ydir','reverse')
        axis equal
        box off

      else
        % Remove flag, we want to draw the grey neurons first
        
        okIdx = find(~isnan(Etot));
        noIdx = find(isnan(Etot));

        set(fig,'currentaxes', axesRetTot);      
        cla, hold on
        plotEnergy(obj.RGCnt(noIdx),obj.RGCdv(noIdx),totCol(noIdx,:),markSize);
        plotEnergy(obj.RGCnt(okIdx),obj.RGCdv(okIdx),totCol(okIdx,:),markSize);
        hold off
        title('Total energy change')
        ylabel('Ventral-Dorsal')
        xlabel('Temporal-Nasal')
        set(gca,'xdir','reverse','ydir','reverse')
        axis equal
        box off
        
        set(fig,'currentaxes', axesRetChem);
        cla, hold on
        plotEnergy(obj.RGCnt(noIdx),obj.RGCdv(noIdx),chemCol(noIdx,:),markSize);
        plotEnergy(obj.RGCnt(okIdx),obj.RGCdv(okIdx),chemCol(okIdx,:),markSize);
        hold off
        title('Chemical energy change')
        ylabel('Ventral-Dorsal')
        xlabel('Temporal-Nasal')
        set(gca,'xdir','reverse','ydir','reverse')
        axis equal
        box off
        
        set(fig,'currentaxes', axesRetAct);    
        cla, hold on
        plotEnergy(obj.RGCnt(noIdx),obj.RGCdv(noIdx),actCol(noIdx,:),markSize);
        plotEnergy(obj.RGCnt(okIdx),obj.RGCdv(okIdx),actCol(okIdx,:),markSize);
        hold off
        title('Activity energy change')
        ylabel('Ventral-Dorsal')
        xlabel('Temporal-Nasal')
        set(gca,'xdir','reverse','ydir','reverse')
        axis equal
        box off
      
        set(fig,'currentaxes', axesRetComp);
        cla, hold on
        plotEnergy(obj.RGCnt(noIdx),obj.RGCdv(noIdx),compCol(noIdx,:),markSize);
        plotEnergy(obj.RGCnt(okIdx),obj.RGCdv(okIdx),compCol(okIdx,:),markSize);
        hold off
        title('Competition energy change')
        ylabel('Ventral-Dorsal')
        xlabel('Temporal-Nasal')
        set(gca,'xdir','reverse','ydir','reverse')
        axis equal
        box off

      end
        
      set(fig,'currentaxes',axesRetTot);
      addColorbar(eMax,eMin,colMap);
      set(fig,'currentaxes',axesRetAct);
      addColorbar(eMax,eMin,colMap);
      set(fig,'currentaxes',axesRetChem);
      addColorbar(eMax,eMin,colMap);
      set(fig,'currentaxes',axesRetComp);
      addColorbar(eMax,eMin,colMap);
 
    else
      % No neuron selected, just clear subplots

      set(fig,'currentaxes', axesRetTot);      
      cla
      
      set(fig,'currentaxes', axesRetChem);
      cla
      
      set(fig,'currentaxes', axesRetAct);    
      cla
      
      set(fig,'currentaxes', axesRetComp);
      cla
            
    end
      
    set(fig,'currentaxes', axesRetMap);
    cla, hold on
    for i = 1:obj.nRGC
      plot(obj.RGCnt(i),obj.RGCdv(i), ...
           '.', 'color', RGCcol(i,:), 'markersize', markSize);
    end
    
    % Mark the synapses
    if(~isempty(selectIdx))
      for i = 1:obj.numPresynapticConnections(selectIdx)
         idx = obj.presynapticConnections(i,selectIdx);
         w = obj.presynapticWeight(i,selectIdx);
         plot(obj.RGCnt(idx),obj.RGCdv(idx), ...
              'yo', 'markersize', 4+w);
      end
    end
    
    hold off
    xlabel('Temporal-Nasal')
    ylabel('Ventral-Dorsal')
    set(gca,'xdir','reverse','ydir','reverse')
        
    title('Retinal map')
    axis equal
    box off

      
    disp(' Done.')
    
  end

  function addColorbar(eMax,eMin,colMap);
  
    cb = colorbar; 
      
    eSpacing = 10^floor(log10(eMax-eMin));
    yValTop = 0:eSpacing:eMax;
    yValBottom = 0:-eSpacing:eMin;
    yVal = unique([yValBottom(end:-1:1), yValTop]);
    
    ytick = interp1(linspace(eMin,eMax,size(colMap,1)), ...
                    1:size(colMap,1),yVal);

    yTickLabel = {};
    for i = 1:numel(yVal)
      yTickLabel{i} = num2str(yVal(i));
    end
    
    set(cb,'ytick',ytick,'yticklabel',yTickLabel) 
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function plotEnergy(RGCdv,RGCnt,col,markSize)
    
    for i = 1:numel(RGCdv)
      plot(RGCdv(i),RGCnt(i),'.', ...
           'color',col(i,:),'markersize',markSize);
    end
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function [idx,dist] = findClosestSC(x,y)
    [dist,idx] = min(sqrt((obj.SCap - x).^2 + (obj.SCml - y).^2));
  end
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function mouseHandle(source, event)

    set(fig,'currentaxes',axesSC);
    xy = get(axesSC,'currentpoint');
    x = xy(1,1);
    y = xy(1,2);

    a = axis();
    
    if(a(1) <= x & x <= a(2) ...
       & a(3) <= y & y <= a(4))

      % Valid point, find closest SC
      [selectIdx,dist] = findClosestSC(x,y);
      
      if(dist > 0.1)
        selectIdx = [];
      end

      plotFigure();
      
    end
    
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
end