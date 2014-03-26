% This function takes a sequence (from loadIterSequence) and
% decimates it, useful when you are playing around with simulations
% that have too many saved data steps.

% If stepsize is a scalar, iterations spaced like that are used, if
% stepsize is a vector, those iterations are extracted
% Note that stepSize is in average iteration per SC neuron, not
% total simulation iterations.

function seqFilt = filterSequence(obj,seq,stepSize,skipZeroFlag,keepFirstNonZeroFlag)
    
  if(isempty(seq))
    seqFilt = [];
    return
  end
    
  availableSteps = cat(1,seq.curStep);

  if(numel(stepSize) > 1)
    requestedSteps = stepSize;
  else
    if(exist('skipZeroFlag') & skipZeroFlag)
      requestedSteps = stepSize:stepSize:max(availableSteps/seq(1).nSC);    
    else
      requestedSteps = 0:stepSize:max(availableSteps/seq(1).nSC);
    end
  end

  requestedSteps = requestedSteps*seq(1).nSC;
  
  if(exist('keepFirstNonZeroFlag') & keepFirstNonZeroFlag)
    nzIdx = find(availableSteps > 0,1,'first');
    requestedSteps = [requestedSteps,availableSteps(nzIdx)];
  end
  
  useStepIdx = zeros(size(requestedSteps));
  
  for i = 1:numel(requestedSteps)
    [~,minIdx] = min(abs(availableSteps - requestedSteps(i)));
    useStepIdx(i) = minIdx(1);
  end
  
  useStepIdx = unique(useStepIdx);
  
  seqFilt = seq(useStepIdx);
  
end