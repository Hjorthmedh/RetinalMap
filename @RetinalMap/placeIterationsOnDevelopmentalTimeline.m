% What iteration corresponds to P0, P1, ... etc
%
% Code currently only handles WT.
%

function placeIterationsOnDevelopmentalTimeline(obj)

 try

  nRep = 200;

  [age,kSigmoidWT,kSigmoidB2KO] = obj.fitSegregationData();

  sequence = obj.loadIterSequence();

  % If the first RetinalMap object in the sequence has no synapses,
  % then throw it out.
  if(numel(sequence) > 0 ...
     & sum(sequence(1).numPresynapticConnections) == 0)
    sequence = sequence(2:end);
    disp('Skipping first object in sequence, it has 0 synapses')
  end

  % Dans developmental data is for WT, AP axis
  % (There is also B2KO and ML axis data, but incomplete)
  
  d = {};
  sep = {};
  kSigmoid = zeros(numel(sequence),1);
  seqIter = cat(1,sequence.curStep)./cat(1,sequence.nSC);
  
  for i = 1:numel(sequence)
    fprintf('%d/%d ',i,numel(sequence))
    [d{i},sep{i},kSigmoid(i)] = ...
        sequence(i).virtualInjectionSegregationExperiment('MLfix',nRep);
  end

  % For each experimental timepoint, find betch matching iteration
  iter = NaN*zeros(size(age));
  iterIdx = NaN*zeros(size(age));  
  kError = NaN*zeros(size(age));
  
  for i = 1:numel(age)
    [matchError,bestMatchIdx] = min(abs(kSigmoid-kSigmoidWT(i)));
    fprintf('age: %d, matchError: %.3f\n', age(i), matchError)
    if(~isempty(bestMatchIdx))
      iterIdx(i) = bestMatchIdx;
      iter(i) = sequence(bestMatchIdx).curStep/sequence(bestMatchIdx).nSC;
      kError(i) = matchError;
    end
  end

  relError = abs(kError./kSigmoidWT);
  
  % Plot the result
  figure
  plot(age,iter,'k-')    
  xlabel('Age (days)')
  ylabel('Iteration')
  title('Dont trust me, unless errors in match checked!')
  
  figure
  subplot(2,1,1)
  plot(age,kSigmoidWT,'r*-')
  xlabel('Age (days)')
  ylabel('kSigmoid')
  
  subplot(2,1,2)
  plot(seqIter,kSigmoid,'k-')
  hold on
  maxK = max(kSigmoid)
  minK = min(kSigmoid);
  a = axis();
  a(4) = min(a(4),10);
  axis(a);
  ofs = (a(4)-a(3))*0.07;
  
  plotIdx = zeros(size(age));
  
  for i = 1:numel(age)
    if(minK <= kSigmoidWT(i) & kSigmoidWT(i) <= maxK)
      plotIdx(i) = 1;
    end
  end

  % Also include the one above and below
  idx = find(kSigmoidWT >= maxK,1,'last');
  if(~isempty(idx))
    plotIdx(idx) = 1;
  end
  
  idx = find(kSigmoidWT <= minK,1,'first');
  if(~isempty(idx))
    plotIdx(idx) = 1;
  end
  
  for i = 1:numel(age)
    if(plotIdx(i))
      xofs = interp1([1 numel(age)],[0.05 0.95],i);
      plot([0 max(seqIter)],[1 1]*kSigmoidWT(i),'r-')
      text(xofs*max(seqIter),kSigmoidWT(i)+ofs, ...
           sprintf('P%d',age(i)), ...
           'color',[1 0 0])
    end
  end
  
  xlabel('Iteration')
  ylabel('kSigmoid')
  
  % keyboard
  
 catch e
   getReport(e)
   keyboard
 end
  
end
