function [xRange,yMedian,yAll] = loessJackKnife(obj,x,y)

  assert(false); % I do not yet trust this function.
                 % Problem: If I send in just the data points, then
                 % I do not get the values of the interpolated
                 % points that I want. If I add them as NaN, then
                 % the span-window is getting confused, because it
                 % counts the NaN points as well for the span.
  
  nGrid = 100;

  span = linspace(20,50,10);
  
  xRange = linspace(0,max(x),nGrid);
  yAll = NaN*zeros(nGrid,numel(x));
  
  for i = 1:numel(x)

    idx = setdiff(1:numel(x),i);
    xUse = [x(idx);transpose(xRange)];
    yUse = [y(idx);transpose(NaN*xRange)]; % The xRange points got no values
                                % to use
    
    yGridL = zeros(numel(xRange),numel(span));
    spanScore = zeros(1,numel(span));
    
    for j = 1:numel(span)
      yL = smooth(xUse,yUse,span(j),'loess');
      
      % Calculate the score from the data point mismatch
      spanScore(j) = sum((yL(1:numel(idx))-yUse(1:numel(idx))).^2);
      
      % Save the y values corresponding to the xRange
      yGridL(:,j) = yL(numel(idx)+1:end);
      
    end

    % Save the y-values corresponding to the best fit, for xRange values
    [~,jBest] = min(spanScore);
    yAll(:,i) = yGridL(:,jBest); 
            
  end
    
  % Just a safety check, there should be no NaNs
  if(1)
    try 
      assert(all(~isnan(yAll(:))));
    catch e
      getReport(e)
      keyboard
    end
  end
    
  yMedian = median(yAll,2);
  
end