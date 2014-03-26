% This function centers all projections and superposes them on top
% of each other

function [H,majorAng,mmRatio,fig] = superposedProjection(obj,RGCidx,showRange)

  if(~exist('RGCidx') | isempty(RGCidx))
    RGCidx = 1:obj.nRGC;
    allRGC = 1;
  else
    allRGC = 0;
  end
  
  if(obj.plotFigures == 0)
    disp('superposedProjection: plotFigures = 0, hiding figures.')
    visFlag = 'off';
  else
    visFlag = 'on';
  end 
  
  if(~exist('showRange'))
    showRange = true;
  end
  
  obj.convertConnectionTables('pre2post');

  [RGCcentAP,RGCcentML] = obj.RGCprojectionCentroids();

  allPoints = [];
  allWeights = [];
  
  for i = 1:numel(RGCidx);
    nCon = obj.numPostsynapticConnections(RGCidx(i));
    idx = obj.postsynapticConnections(1:nCon,RGCidx(i));
    
    ap = obj.SCap(idx) - RGCcentAP(RGCidx(i));
    ml = obj.SCml(idx) - RGCcentML(RGCidx(i));
    
    allPoints = [allPoints; ap, ml];
    allWeights = [allWeights; obj.postsynapticWeight(1:nCon,i)];
    
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if(0)
    figure
    hold on
    for i = 1:size(allPoints,1)
      plot(allPoints(i,1),allPoints(i,2),'.','markersize',3+allWeights(i));
    end
    hold off
    xlabel('Anterior - Posterior')
    ylabel('Medial - Lateral')
  end

  clear ap ml
  
  if(showRange)
    % Only do ellipses for RGC with synapses!
    idxOk = find(obj.numPostsynapticConnections(RGCidx) > 0);
    
    % Make post stamp plot to verify range of variations
    for i = 1:numel(idxOk)
      nCon = obj.numPostsynapticConnections(idxOk(i));
      idx = obj.postsynapticConnections(1:nCon,idxOk(i));
      
      ap{i} = obj.SCap(idx) - RGCcentAP(idxOk(i));
      ml{i} = obj.SCml(idx) - RGCcentML(idxOk(i));
      w = obj.postsynapticWeight(1:nCon,idxOk(i));
      
      [xeS{i},yeS{i},mmRatioS(i),majorAngS(i)] = makeEllipse(ap{i}, ml{i}, w, 100);
    end
    
    [~,idx] = sort(mmRatioS);
    
    showIdx = idx(round(linspace(1,numel(idxOk),16)));
    
    fig(2) = figure('visible',visFlag);

    for i = 1:numel(showIdx)
      sp(i) = subplot(4,4,i);

      plot(ap{showIdx(i)},ml{showIdx(i)},'o','color',[1 1 1]*0.6);
      hold on
      plot(xeS{showIdx(i)},yeS{showIdx(i)},'k-');
      hold off
      box off
      axis equal
      axis off
      title(sprintf('%.3f', mmRatioS(showIdx(i))));
    end
    
    % Set all the axis equal
    aMin = [inf inf];
    aMax = [-inf -inf];
    
    for i = 1:numel(sp)
      set(fig(2),'currentaxes',sp(i));
      a = axis();
      aMin = min(a([1 3]),aMin);
      aMax = max(a([2 4]),aMax);
    end

    aNew = [aMin(1) aMax(1) aMin(2) aMax(2)];
    
    for i = 1:numel(sp)
      set(fig(2),'currentaxes',sp(i));
      axis(aNew);
    end
    
    
    
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  binRange = -0.25:0.01:0.25;
  nBins = numel(binRange);
  
  X = round(interp1(binRange,1:nBins,allPoints(:,1),'linear','extrap'));
  Y = round(interp1(binRange,1:nBins,allPoints(:,2),'linear','extrap')); 
  
  okIdx = find(1 <= X & X <= nBins & 1 <= Y & Y <= nBins);
  
  % 2D histogram
  H = accumarray([Y(okIdx), X(okIdx)],allWeights(okIdx),[nBins nBins]);
  
  if(0)
    % Just double checking accumarray
    H2 = zeros(nBins,nBins);
    for i = 1:numel(okIdx)
      ik = okIdx(i);
      H2(Y(ik),X(ik)) = H2(Y(ik),X(ik)) + allWeights(ik);
    end
  end
    

  
  fig(1) = figure('visible',visFlag);
  imagesc(binRange,binRange,H)
  hold on
  nPoints = 100;
  [xe,ye,mmRatio,majorAng] = makeEllipse(allPoints(:,1), ...
                                         allPoints(:,2), ...
                                         allWeights, ...
                                         nPoints);

  plot(xe,ye,'b-','linewidth',2)
  
  hold off
  box off
  
  % title(sprintf('%s %d synapses', obj.simName, sum(H(:))))
  if(allRGC)
    title(sprintf(['Superposed projections %s,\n' ...
                   'major angle %.3f, ratio %.3f'], ...
                  obj.simName, majorAng*180/pi, mmRatio),'fontsize',24)  
  else
    title(sprintf(['Superposed projections (sub population: %d),\n' ...
                   'major angle %.3f, ratio %.3f'], ...
                  numel(RGCidx), ...
                  majorAng*180/pi, mmRatio),'fontsize',24)  
    
  end
  xlabel('Anterior - Posterior','fontsize',20)
  ylabel('Medial - Lateral','fontsize',20)
  set(gca,'fontsize',16)
  colormap('hot')
  axis equal
  
  figName = sprintf('%s/%s-superposed-projection.pdf', ...
                    obj.figurePath, obj.simName);
  
  switch(obj.plotFigures)
    case {0,1}
      saveas(fig(1),figName,'pdf');
    case 2
      obj.addParametersToFigure(strrep(figName,'.pdf','.eps'));
  end
 
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function [covMat,mu] = weightedCov(x,y,w)

    xMean = sum(x.*w)/sum(w);
    yMean = sum(y.*w)/sum(w);
    
    covMat(1,1) = sum((x-xMean).*w .* (x-xMean).*w)/sum(w.*w);
    covMat(1,2) = sum((x-xMean).*w .* (y-yMean).*w)/sum(w.*w);
    covMat(2,1) = covMat(1,2);
    covMat(2,2) = sum((y-yMean).*w .* (y-yMean).*w)/sum(w.*w);
    
    mu = [xMean,yMean];
    
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function [xe,ye,mmRatio,majorAng] = makeEllipse(ap,ml,w,nPoints)
  
    [covMat,mu] = weightedCov(ap,ml,w);
  
    [V,D] = eig(covMat);
      
    ang = linspace(0,2*pi,nPoints);
    xc = cos(ang);
    yc = sin(ang);
    
    xy = [transpose(xc),transpose(yc)]*sqrt(D)*V;
    
    xe = xy(:,1) + mu(1);
    ye = xy(:,2) + mu(2);
    
    [~,majorIdx] = max(diag(D));
    [~,minorIdx] = min(diag(D));
    
    % Major and minor axis
    mmRatio = sqrt(D(majorIdx,majorIdx)/D(minorIdx,minorIdx));
    majorAng = atan(V(2,majorIdx)/V(1,majorIdx));
    
  end
    
end