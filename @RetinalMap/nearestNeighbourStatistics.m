function [nearestNeighbourRegularityIndex, voronoiRegularityIndex] = ...
      nearestNeighbourStatistics(obj,structure,X,Y)

  switch(structure)
    case {'Retina','retina'}
      d = obj.retinalDistance(); 
    case {'SC','sc'}
      d = obj.SCDistance();
    otherwise
      fprintf('nearestNeighbourStatistics: Unknown structure %s\n', structure)
      keyboard
  end
  
  try
    d = d + diag(max(d(:))*ones(size(d,1),1));
  catch e
    getReport(e)
    keyboard
  end
    
  dMin = min(d);
  
  minDmean = mean(dMin);
  minDstd = std(dMin);
  
  nearestNeighbourRegularityIndex = minDmean/minDstd;
  
  if(0)
  
    figure
    hist(dMin(:),50);
    hold on
    a = axis();
    plot(minDmean*[1 1],a(3:4),'r-', ...
         (minDmean+minDstd)*[1 1],a(3:4),'r--', ...
         (minDmean-minDstd)*[1 1],a(3:4),'r--');
    hold off
    xlabel('Closest neighbour distance')
    ylabel('Count')

    title(sprintf('Nearest neighbour regularity index %.2f', ...
                  nearestNeighbourRegularityIndex))
  
  end
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    switch(structure)
      case {'Retina','retina'}
        X = [obj.RGCnt,obj.RGCdv];
      case {'SC','sc'}
        X = [obj.SCap, obj.SCml];
      otherwise
        fprintf('nearestNeighbourStatistics: Unknown structure %s\n', structure)
        keyboard
    end

    % Calculate the voronoi regularity index

    % We want to exclude cells that are outside
    K = convhull(X(:,1),X(:,2));
    Kx = X(K,1);
    Ky = X(K,2);
    
    [V,C] = voronoin(X);
    
    vorArea = NaN*zeros(numel(C),1);
    for i = 1:numel(C)
      x = V(C{i},1);
      y = V(C{i},2);

      if(all(inpolygon(x,y,Kx,Ky)))
        vorArea(i) = polyarea(x,y);
      end
    end
    
    idx = find(~isnan(vorArea));
    
    vorArea = vorArea(idx);

    voronoiRegularityIndex = mean(vorArea)/std(vorArea);
    
    if(0)
      figure
      hist(vorArea,50)
      xlabel('Voronoi area')
      ylabel('Count')
      title(sprintf('Voronoi area regularity index %.2f', ...
                    voronoiRegularityIndex))
    end
    
end
