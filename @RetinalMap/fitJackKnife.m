function [k,kAll,xRange,yMedian,yAll] = fitJackKnife(obj,x,y,func,k0)

  if(~exist('k0'))
    k0 = [];
  end
  
  nGrid = 100;
  kAll = zeros(2,numel(x));

  xRange = linspace(0,max(x),nGrid);
  yAll = zeros(nGrid,numel(x));
  
  k = obj.fitSegregation(x,y,k0);
  
  for i = 1:numel(x)

    idx = setdiff(1:numel(x),i);
    xUse = x(idx);
    yUse = y(idx);
    
    kAll(:,i) = obj.fitSegregation(xUse,yUse,k0);

    yAll(:,i) = func(kAll(:,i),xRange);
    
  end
    
  % Just a safety check, there should be no NaNs
  try 
    assert(all(~isnan(yAll(:))));
  catch e
    getReport(e)
    keyboard
  end
    
  yMedian = median(yAll,2);
  
end