% yLow is 2.5 percentile

% yHigh is 97.5 percentile
function [xRange,yMedian,yLow,yHigh,yAll] = fitErrorEstimateN(obj,x,y,func,k0,nRep,nPoints)
    
  if(~exist('nRep'))
    nRep = 1000;
  end
  
  % Number of points for the curve
  if(~exist('nPoints'))
    nPoints = 100;
  end
  
  xRange = transpose(linspace(min(x),max(x),nPoints));
  
  idxLow = ceil(0.025*nRep);
  idxHigh = ceil(0.975*nRep);
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Original fit
  k = nlinfit(x,y,func,k0);
  
  kAll = zeros(numel(k),nRep);
  yAll = zeros(nPoints,nRep);
  
  for iRep = 1:nRep
    
    % 1. Resample data
    [xNew,yNew] = resampleData(x,y);
    
    % 2. Fit function to the distribution
    kAll(:,iRep) = nlinfit(xNew,yNew,func,k0);
    yAll(:,iRep) = func(kAll(:,iRep),xRange);
    
  end

  % Mean fit
  yMedian = median(yAll,2);

  % Find the parameter range for 95% coverage 
  yAllSorted = sort(yAll,2);
  
  yLow = yAllSorted(:,idxLow);
  yHigh = yAllSorted(:,idxHigh);

  %keyboard  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function [xNew,yNew] = resampleData(x,y)

    idx = ceil(numel(x)*rand(size(x,1),size(x,2)));
    xNew = x(idx);
    yNew = y(idx);
  
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
end