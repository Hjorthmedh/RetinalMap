function fig = plotTerminationZoneSize(obj,sequence,percentile,titleText)

  skipTables = true;
  
  if(~exist('sequence') | isempty(sequence))
    sequence = obj.loadIterSequence(skipTables); 
  end
 
  if(~exist('percentile'))
    percentile = 0.95;
  end
  
  if(numel(sequence) > 0 ...
     & sum(sequence(1).numPresynapticConnections) == 0)
    sequence = sequence(2:end);
    disp('Skipping first object in sequence, it has 0 synapses')
  end
  
  if(isempty(sequence))
    disp('No sequence loaded, or no synapses found.')
    fig = [];
    return
  end
 
  iter = zeros(numel(sequence),1);
  spread95Median = zeros(numel(sequence),numel(percentile));
  spread95Per5 = zeros(numel(sequence),numel(percentile));
  spread95Per95 = zeros(numel(sequence),numel(percentile));      

  nRGCLabeledMean = zeros(numel(sequence),1);
  nRGCLabeledSEM  = zeros(numel(sequence),1);  
  nSCLabeledMean  = zeros(numel(sequence),1);
  nSCLabeledSEM   = zeros(numel(sequence),1);  
  
  for i = 1:numel(sequence)   
    iter(i) = sequence(i).curStep/sequence(i).nSC;
    fprintf('Analysing file %d/%d (iter %d)\n', i, numel(sequence), iter(i))
    % spread = sequence(i).calculateTerminationZoneSize();
    [spread,labeledRGC,labeledSC] = ...
        sequence(i).calculateTerminationZoneSizeInjection([],percentile);    

    spread95Median(i,:)  = median(spread);
    spread95Per5(i,:)    = prctile(spread,5);
    spread95Per95(i,:)   = prctile(spread,95);
    
    nRGCLabeledMean(i) = mean(labeledRGC);
    nRGCLabeledSEM(i)  = std(labeledRGC)/sqrt(numel(labeledRGC)-1);    
    nSCLabeledMean(i)  = mean(labeledSC);
    nSCLabeledSEM(i)   = std(labeledSC)/sqrt(numel(labeledSC)-1);
    
  end

  if(obj.plotFigures == 0)
    disp('plotTerminationZoneSize: plotFigures = 0, hiding figures!')
    visFlag = 'off';
  else
    visFlag = 'on';
  end

  
try  
  
  % Plot the termination zone size
  fig(1) = figure('visible',visFlag);
  
  for i = 1:numel(percentile)
    idx = find(~isnan(spread95Per5(:,i)));
    a = area(repmat(iter(idx),1,2), ...
             [spread95Per5(idx,i), ...
              spread95Per95(idx,i)-spread95Per5(idx,i)], ...
             'facecolor', 0.6*[1 1 1], ...
             'edgecolor', 0.6*[1 1 1]);
  
    delete(a(1));
    alpha(get(a(2),'children'),0.5);

    hold all
      
    % Switch to mean + SEM instead?
    p(i) = plot(iter,spread95Median(:,i),'-','linewidth', 5);
    
    pLeg{i} = sprintf('Percentile %d', 100*percentile(i))
    
    xlabel('Iterations','fontsize',20)
    % ylabel('Retinal spread (radius)','fontsize',20)   % !!! old
    ylabel('Fraction of retina marked','fontsize',20)    
    set(gca,'fontsize',18)
    box off
  end
  legend(p,pLeg)
  
  if(exist('titleText'))
    title(titleText,'fontsize',25)
  end
  
  fName = sprintf('%s/%s-termination-zone-spread.pdf', ...
                  sequence(1).figurePath, sequence(1).simName);

  switch(obj.plotFigures)
    case {0,1}
      saveas(fig(1),fName,'pdf');
    case 2
      obj.addParametersToFigure(strrep(fName,'.pdf','.eps'));
  end
  
  
  % Plot the number of RGC neurons labeled
  fig(2) = figure('visible',visFlag);
  
  idx = find(~isnan(nRGCLabeledMean));
  a = area(repmat(iter(idx),1,2), ...
           [nRGCLabeledMean(idx) - nRGCLabeledSEM(idx), ...
            2*nRGCLabeledSEM(idx)], ...
            'facecolor', 0.6*[1 1 1], ...
            'edgecolor', 0.6*[1 1 1]);

  delete(a(1));
  hold on
  
  p = plot(iter,nRGCLabeledMean,'k-','linewidth',5);
  
  xlabel('Iterations','fontsize',20)
  ylabel('#RGC labeled','fontsize',20)    
  set(gca,'fontsize',18)
  box off
    
  title(sprintf('%s: showing mean and SEM', obj.simName))
    
  
  fName = sprintf('%s/%s-number-of-RGC-labeled.pdf', ...
                  sequence(1).figurePath, sequence(1).simName);

  switch(obj.plotFigures)
    case {0,1}
      saveas(fig(2),fName,'pdf');
    case 2
      obj.addParametersToFigure(strrep(fName,'.pdf','.eps'));
  end  
  
catch e
    getReport(e)
    keyboard
end
  
      
end