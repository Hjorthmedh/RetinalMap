% Find the iteration in the sequence for the current object that
% best matches matchObj segregation curve.

% If kMatch is given, matchObj is ignored 
% If kMatch is not given, it is calculated from matchObj

% The option to pass kMatch exists so that all can be compared to the
% same reference curve.

function [bestIter,bestIdx,kFit,kMatch] = findSegregationCrossing(obj,matchObj,kMatch,stopAtIter)

  debugPlots = true;
  
  skipTable = 2; % Conserve more memory
  seq = obj.loadIterSequence(skipTable);
  
  if(sum(seq(1).totalWeightRGC) == 0)
    % Discard t=0 if it has no synapses
    seq = seq(2:end);
  end

  if(exist('stopAtIter') & ~isempty(stopAtIter))
    keepIdx = find(cat(1,seq.curStep) <= stopAtIter*seq(1).nSC);
    seq = seq(keepIdx);
  end
  
  nInj = 200;
  injectionType = 'MLfix';
  smoothWindow = 10;

  fitScore = zeros(numel(seq),1);  
  kFit = zeros(2,numel(seq));
  
  % Calculate all the fitting curves, nSets, each with nInj injections
  for i = 1:numel(seq)
    fprintf('Seq: %d/%d\n', i, numel(seq))
    
    [distance,segregation,kFit(:,i)] = ...
        seq(i).virtualInjectionSegregationExperiment(injectionType,nInj);
  end
  
  xRange = linspace(0,1,50);
  
  if(~exist('kMatch') | numel(kMatch) ~= 2)
    [~,~,kMatch] = matchObj.virtualInjectionSegregationExperiment(injectionType,nInj);    
  end
  
  yMatch = logistic(kMatch,xRange);
  
 try
  
  distTotal = zeros(numel(seq),1);
  crossFlag = false;  
  
  % We will only report the closest match if the curve is able to
  % cross the reference line at some point. If it never reaches it,
  % then we report it as NaN
  
  for i = 1:numel(seq)
        
    y = logistic(kFit(:,i),xRange);
    distTotal(i) = sum((yMatch-y).^2);
    
    if(sum(y-yMatch) > 0)
      crossFlag = true;
    end
    
  end

  distTotalSmoothed = smooth(distTotal,smoothWindow);
  
  [~,idx] = min(distTotalSmoothed);

  bestIdx = idx(1);      
  bestIter = seq(idx(1)).curStep / seq(idx(1)).nSC;      

  if(debugPlots)
    figure
    subplot(2,1,1)
    hold on
        
    % Draw the fit line
    y = logistic(kFit(:,bestIdx),xRange);
    yRef = logistic(kMatch,xRange);
    p = plot(xRange,y,'k-',xRange,yRef,'--');
    legend(p,'Best match','Reference')
    set(gca,'fontsize',15)
    xlabel('Normlised distance on SC','fontsize',20)
    ylabel('Retinal segregation','fontsize',20)
    box off
    
    if(crossFlag)
      title('Ok fit')
    else
      title('Fails to reach segregation')
    end
    
    subplot(2,1,2)
    iter = cat(1,seq.curStep) ./ cat(1,seq.nSC);
    plot(iter,distTotal,'k-', ...
         iter,distTotalSmoothed,'r-', ...
         iter(bestIdx),distTotalSmoothed(bestIdx),'r*');
    xlabel('Iteration','fontsize',20)
    ylabel('Mismatch','fontsize',20)
    box off
    set(gca,'fontsize',15)
    
  end
  
  
  if(~crossFlag)
    bestIter = NaN;
    bestIdx = NaN;
  end
    
    
 catch e
   getReport(e)
   keyboard
 end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function y = logistic(k,x)
    % y = -0.5./(1+(x/inflection).^slope)+1;
    y = 1 - 0.5./(1+(x/k(1)).^k(2));
    
    % If k(1) is negative, we can get imaginary solutions
    % If k(2) is negative, x=0 --> y = 1.
    if(any(k < 0))
      y = ones(size(x))*1e5;
    end
    
  end
  
  
end