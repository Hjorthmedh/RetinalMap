function inspectSelfAddCost(obj)

  % Find all the RGC centroids for the SC neurons

  centroidRGCidx = zeros(obj.nSC,1); 
      
  for j = 1:obj.nSC
    centroidRGCidx(j) = getRGCcentroidNeuron(j); 
  end
      
  % For each SC neuron, what would be the cost of adding a synapse to that RGC

  Etot = zeros(obj.nSC,1);
  Echem = zeros(obj.nSC,1);
  Eact = zeros(obj.nSC,1);
  Ecomp = zeros(obj.nSC,1);
            
  for i = 1:obj.nSC
    [Etot(i),Echem(i),Eact(i),Ecomp(i)] = ...
        obj.calculateSynapseAdditionEnergy(i,centroidRGCidx(i));
  end
      
  % Set up a colormap    
  eMax = max([max(Etot),max(Echem),max(Eact),max(Ecomp)]);
  eMin = min([min(Etot),min(Echem),min(Eact),min(Ecomp)]);
      
  colormap('jet');
  colMap = colormap('jet');
      
  eSpace = linspace(eMin,eMax,64);
  
  totCol = zeros(obj.nSC,3);
  chemCol = zeros(obj.nSC,3);
  actCol = zeros(obj.nSC,3);
  compCol = zeros(obj.nSC,3);   
  
  for i = 1:3
    totCol(:,i) = interp1(eSpace,colMap(:,i),Etot);
    chemCol(:,i) = interp1(eSpace,colMap(:,i),Echem);
    actCol(:,i) = interp1(eSpace,colMap(:,i),Eact);
    compCol(:,i) = interp1(eSpace,colMap(:,i),Ecomp);     
  end
  
  % Plot the resulting cost
  figure, hold on
  for i = 1:obj.nSC
    plot(obj.SCml(i),obj.SCap(i),'.', 'color', totCol(i,:), ...
         'markersize', 12);
  end
  xlabel('Medio-Lateral','fontsize',24)
  ylabel('Anterior-Posterior','fontsize',24)
  title('Energy to add synapse to RGC centroid',...
        'fontsize',24)
  set(gca,'fontsize',20)
  addColorbar(eMax,eMin,colMap);

  saveas(gcf,sprintf('%s/%s-cost-add-synapse-to-centroid-iter-%d.pdf', ...
                     obj.figurePath, obj.simName, obj.curStep/obj.nSC),'pdf')
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  function RGCidx = getRGCcentroidNeuron(SCidx)
  
    dvC = 0; ntC = 0;
      
    % Find the coordinates of the centroid
    for i = 1:obj.numPresynapticConnections(SCidx)
      idx = obj.presynapticConnections(i,i);
      w = obj.presynapticWeight(i,i);
      dvC = dvC + w*obj.RGCdv(idx);
      ntC = ntC + w*obj.RGCnt(idx);
    end

    if(obj.numPresynapticConnections(SCidx) > 0)
      dvC = dvC / obj.totalWeightSC(SCidx);
      ntC = ntC / obj.totalWeightSC(SCidx);
    else
      dvC = NaN;
      ntC = NaN;
    end
      
    % Find the closest matching RGC
    [~,RGCidx] = min(sqrt((obj.RGCnt-ntC).^2 + (obj.RGCdv-dvC).^2));
    
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
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

  
end
