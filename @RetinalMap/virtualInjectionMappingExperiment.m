% This experiment does a set of injections, to determine the
% mapping from retina to SC, NT to AP, or DV to ML.
%
% Retina -> SC

function fig = virtualInjectionMappingExperiment(obj,axisDir)

  nInj = 10;
  
  x = linspace(0,1,nInj) + 0.1*rand(1,nInj);
  x = min(max(x,0),1);
  
  switch(axisDir)
    case {'AP','ap','NT','nt'}
      NT = x*(max(obj.RGCnt)-min(obj.RGCnt)) + min(obj.RGCnt);
      DV = 0.5*ones(1,nInj);
        
      xLabel = 'Nasal - Temporal';
      yLabel = 'Anterior - Posterior';
      
    case {'ML','ml','DV','dv'}
      NT = 0.5*ones(1,nInj);
      DV = x*(max(obj.RGCdv)-min(obj.RGCdv)) + min(obj.RGCdv);
      
      xLabel = 'Dorsal - Ventral';
      yLabel = 'Medial - Lateral';
      
    otherwise
      fprintf('Unknown axis: %s\n', axisDir)
      keyboard
  end

  % Check if injection site is ok
  retinaWidth = max(obj.RGCnt)-min(obj.RGCnt);
  injRadius = retinaWidth/2 * sqrt(0.01); % Retinal injection size
  
  % Reber (email 10 oct 2012) estimate 1-3% of retinal area
  % marked. If we assume that the retina is a disk, 
  % 0.01*pi*R^2 = pi*r^2 (r = inj radius, R = retinal radius)
  % r = R*sqrt(0.01), ie r = 0.05

  
    
  X = [];
  Y = [];
  
  Xcent = [];
  Ycent = [];
  
  for i = 1:numel(NT)
    
    % Find RGC close to injection site
    preIdx = find((obj.RGCnt-NT(i)).^2 + (obj.RGCdv-DV(i)).^2 < injRadius^2);

    % Find what SC neurons they project to
    postIdx = findPostsynapticSC(preIdx);
    
    if(isempty(postIdx))
      % This injection had no neurons close to it, skipping
      fprintf('Found no neurons with synapses close to injection site (%.2f,%,2f)\n', ...
              NT(i),DV(i))
      continue
    end
    
    switch(axisDir)
      case {'AP','ap','NT','nt'}
        [nt,ap] = meshgrid(obj.RGCnt(preIdx),obj.SCap(postIdx));
          
        X = [X; nt(:)];
        Y = [Y; ap(:)];
        
        pc = findProjectionCentres(obj.SCap(postIdx));
        Xcent = [Xcent; mean(obj.RGCnt(preIdx))*ones(size(pc))];
        Ycent = [Ycent; pc];
        

      case {'ML','ml','DV','dv'}
        [dv,ml] = meshgrid(obj.RGCdv(preIdx),obj.SCml(postIdx));
        
        X = [X; dv(:)];
        Y = [Y; ml(:)];

        Xcent = [Xcent; mean(obj.RGCdv(preIdx))*ones(size(pc))];
        Ycent = [Ycent; pc];
        
    end

  end

  % Make the plots
  
  if(obj.plotFigures == 0)
    disp('virtualInjectionMappingExperiment: plotFigures = 0, hiding figures!')
    visFlag = 'off';
  else
    visFlag = 'on';
  end  
  
  fig = figure('visible',visFlag);

  plot(X,Y,'.','color',[1 1 1]*0.6)
  hold on
  plot([0 1],[1 0],'k--'), hold on
  plot(Xcent,Ycent,'k.','markersize',15)
  xlabel(xLabel,'fontsize',20)
  ylabel(yLabel,'fontsize',20)
  set(gca,'fontsize',20)
  box off
  axis equal
  axis([0 1 0 1])
  
  fName = sprintf('%s/%s-virtual-exp-%s-axis.pdf', ...
                  obj.figurePath, obj.simName, axisDir);
  
  switch(obj.plotFigures)
    case {0,1}
      saveas(gcf,fName,'pdf');
    case 2
      obj.addParametersToFigure(strrep(fName,'.pdf','.eps'));
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function postIdx = findPostsynapticSC(RGCidx)
    
    if(isempty(RGCidx))
      postIdx = [];
      return
    end
    
    % Make our life easier...
    obj.convertConnectionTables('pre2post');
  
    nMax = sum(obj.numPostsynapticConnections(RGCidx));
    nextIdx = 1;
    
    for j = 1:numel(RGCidx)
      nCon = obj.numPostsynapticConnections(RGCidx(j));
      postIdx(nextIdx:(nextIdx+nCon-1)) = ...
        obj.postsynapticConnections(1:nCon,RGCidx(j));
      
      nextIdx = nextIdx + nCon;
    end
    
    postIdx = unique(postIdx);
    
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function projCenters = findProjectionCentres(SCcoord)
    
    % Lets use k-means to determine if they are ectopic or not, same
    % method as in locateCollapsePoint.m
    
    idx = kmeans(SCcoord,2,'emptyaction','singleton');

    idx1 = find(idx == 1);
    idx2 = find(idx == 2);
    
    n1 = numel(idx1);
    n2 = numel(idx2);
    
    nRatio = min(n1,n2)/max(n2,n2);
    
    mean1 = mean(SCcoord(idx1));
    std1  = std(SCcoord(idx1));
    mean2 = mean(SCcoord(idx2));
    std2  = std(SCcoord(idx2));

    separationReq = 2; %1.5; 
    minRatio = 0.05; % This is seems to be pretty permissive
   
    if(abs(mean1-mean2) > (std1+std2)*separationReq & nRatio > minRatio) 
      projCenters = [mean1; mean2];
    else
      projCenters = mean(SCcoord);
    end    
    
  end
  
end