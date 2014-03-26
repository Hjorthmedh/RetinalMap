% sequences is a cell array of sequences, legends are the
% corresponding legends to use

function fig = plotTerminationZoneSizeCompare(obj,sequences,percentile,legendText,colours)

  skipTables = true;
  fig = [];
  iter = {};
  
  if(~exist('sequences'))
    sequences{1} = obj.loadIterSequence(skipTables); 
  end

  if(~exist('percentile'))
    percentile = 0.95;
  end

  nRep = 50;
  
  for is = 1:numel(sequences)
    if(numel(sequences{is}) > 0 ...
       & sum(sequences{is}(1).numPresynapticConnections) == 0)
      sequence = sequences{is}(2:end);
      disp('Skipping first object in sequence, it has 0 synapses')
    else
      sequence = sequences{is};
    end    
    
    if(isempty(sequence))
      disp('No sequence loaded, or no synapses found.')
      continue
    end
    
    iter{is} = zeros(numel(sequence),1);
    spread95Median{is} = zeros(numel(sequence),1);
    spread95Per5{is} = zeros(numel(sequence),1);
    spread95Per95{is} = zeros(numel(sequence),1);      

    for i = 1:numel(sequence)   
      iter{is}(i) = sequence(i).curStep/sequence(i).nSC;
      fprintf('Analysing file %d/%d (iter %d)\n', i, numel(sequence), iter{is}(i))
      spread = sequence(i).calculateTerminationZoneSizeInjection(nRep,percentile);    
      spread95Median{is}(i) = median(spread);
      spread95Per5{is}(i) = prctile(spread,5);
      spread95Per95{is}(i) = prctile(spread,95);            
    end    
    
  end

  if(isempty(iter))
    % If this is empty, we have not loaded any sequences to display
    return
  end
  
  if(obj.plotFigures == 0)
    disp('plotTerminationZoneSize: plotFigures = 0, hiding figures!')
    visFlag = 'off';
  else
    visFlag = 'on';
  end

  fig = figure('visible',visFlag);
    
  for is = 1:numel(iter)
    % Plot the spread first
    idx = find(~isnan(spread95Per5{is}));
    a = area(repmat(iter{is}(idx),1,2), ...
             [spread95Per5{is}(idx), ...
              spread95Per95{is}(idx)-spread95Per5{is}(idx)], ...
             'facecolor', colours(is,:)*0.4+0.6*[1 1 1], ...
             'edgecolor', colours(is,:)*0.4+0.6*[1 1 1]);
  
    delete(a(1)); 
    hold on
  end
  
  for is = 1:numel(iter)
    p(is) = plot(iter{is},spread95Median{is},'-','linewidth', 5,'color',colours(is,:));
    xlabel('Iterations','fontsize',20)
    % ylabel('Termination zone spread (radius)','fontsize',20)
    ylabel('Projection zone area (fraction)','fontsize',20)  
    set(gca,'fontsize',18)
    box off
  end  

  legend(p,legendText)
  
end