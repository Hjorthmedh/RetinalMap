% This function uses kernel density estimate to analyse the contours
% areaCoverage is in fraction of the retina

% Test:
% r.contourAnalysisKDE(r.makeVirtualInjection(0.5,0.5),[0.95 0.75 0.5 0.25 0.05],true)

function [areaCoverage,kEst,fig] = contourAnalysisKDE(obj, RGCidx, ...
                                                      percentile, plotFig)

  if(~exist('plotFig'))
    plotFig = false;
  end

  if(isempty(RGCidx))
    disp('No RGC cells were selected')
    areaCoverage = zeros(size(percentile));
    kEst = NaN;
    fig = [];
    return
  end
  
  nGrid = 100;
  ntRange = linspace(0,1,nGrid);
  dvRange = linspace(0,1,nGrid);
  [NT,DV] = meshgrid(ntRange,dvRange);
  
  % First, we need to estimate the kernel size. This is done using
  % leave-one-out.
  
  d2 = (kron(obj.RGCnt(RGCidx),ones(1,numel(RGCidx))) ...
        - kron(ones(numel(RGCidx),1),transpose(obj.RGCnt(RGCidx)))).^2 ...
       + (kron(obj.RGCdv(RGCidx),ones(1,numel(RGCidx))) ...
          - kron(ones(numel(RGCidx),1),transpose(obj.RGCdv(RGCidx)))).^2;

  % Set diagonal elements to inf, to avoid them giving a
  % contribution to cross-validation.
  d2 = d2 + diag(inf*(1+diag(d2)));
  
  kEst = fminsearch(@kernelMismatch,0.1);
  
  density = KDEst(kEst,NT,DV);
  
  % Remove the density bits outside the retina
  density = maskRetinaF(density);
  
  % Normalise
  % Should be close to 10000 for a 100x100 grid
  % fprintf('Debug density sum: %f\n', sum(density(:)))
  density = density / sum(density(:));
  
  maxArea = nnz(maskRetinaF());

  % How many points on the grid do we need to get percentile area
  % coverage? We start with the big contributions, and keep adding
  % the largest available one until we reach percentile.
  areaCoverage = zeros(size(percentile));
  densityThreshold = zeros(size(percentile)); % Used for plotting contours
  
  [sortedDensity,sortIdx] = sort(density(:),'descend');
  cumDensity = cumsum(sortedDensity);
  
  for i = 1:numel(percentile)
    idx = find(cumDensity >= percentile(i),1,'first');
    
    if(isempty(idx))
      % When we do masking to only look inside the retina a bit of
      % the probability mass outside is removed, this could mean
      % that we can not get up to the percentile level. Then say
      % that we need the entire retina.
      disp('This should not happen!')
      keyboard
      areaCoverage(i) = 1;
      densityThreshold = 1e-20;
    else
      areaCoverage(i) = idx/maxArea;
      densityThreshold(i) = density(sortIdx(idx));
    end

  end
  
  if(plotFig)
    col = [168, 221, 181; 
           123, 204, 196; 
           78, 179, 211; 
           43, 140, 190; 
           8, 88, 158]/255;

    fig = figure();
    imagesc([0 1],[0 1],density);
    colorbar
    hold on

    for i = 1:numel(percentile)
      [c,h] = contour(NT,DV,density, densityThreshold(i));
      child = get(h,'children');
      set(child,'facecolor','none','edgecolor',col(i,:));
      set(h,'linecolor',col(i,:))
      
      p(i) = h;
      pLeg{i} = sprintf('%d%%',percentile(i)*100);
    end
    legend(p,pLeg);
  
    plot(obj.RGCnt(RGCidx),obj.RGCdv(RGCidx),'w.')
 
    
    xlabel('Temporal - Nasal')
    ylabel('Ventral - Dorsal')
    set(gca, 'xdir','reverse','ydir','reverse');
    axis equal
    axis([0 max(max(obj.RGCnt),1) 0 max(max(obj.RGCdv),1)])
    box off
    
  else
    fig = [];
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function e = kernelMismatch(kEst)
    % Normally we would maximize the cross-validated
    % log-likelyhood, but we are using fminsearch so I put a minus
    % sign in front, so we can minimize.
      
    f = -1/(2*kEst^2);      
    c = 1/(2*pi*kEst^2) / (numel(RGCidx)-1); % -1 from cross-validation

    e = -sum(log(c*sum(exp(f*d2),1)));
    
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function v = KDEst(kEst,X,Y)

    % d2 is distance between RGC (RGCidx) and the grid points (X,Y)
    % RGC on rows, grid on columns
    
    d2 = (kron(obj.RGCnt(RGCidx),ones(1,numel(X))) ...
         - kron(ones(numel(RGCidx),1),transpose(X(:)))).^2 ...
         + (kron(obj.RGCdv(RGCidx),ones(1,numel(Y))) ...
            - kron(ones(numel(RGCidx),1),transpose(Y(:)))).^2;
         
    f = -1/(2*kEst^2);
    c = (1/(2*pi*kEst^2)) / numel(RGCidx);
    
    v = c*sum(exp(f*d2),1);
    
    v = reshape(v,size(X,1),size(X,2));
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function Fmask = maskRetinaF(F)
    
    cIdx = convhull(obj.RGCnt,obj.RGCdv);
    
    [x,y] = convCoords(obj.RGCnt(cIdx),obj.RGCdv(cIdx));

    if(exist('F') & ~isempty(F))
      Fmask = F.*double(poly2mask(x,y,size(F,2),size(F,1)));
    else
      Fmask = poly2mask(x,y,nGrid,nGrid);
    end
    
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  function [x,y] = convCoords(nt,dv)
    
    padding = 0;
    
    x = nGrid*(nt-min(obj.RGCnt(:))+padding)/(max(obj.RGCnt(:))-min(obj.RGCnt(:))+padding*2);
    y = nGrid*(dv-min(obj.RGCdv(:))+padding)/(max(obj.RGCdv(:))-min(obj.RGCnt(:))+padding*2);
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
end
