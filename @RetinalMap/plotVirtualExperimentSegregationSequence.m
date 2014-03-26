function [k,iter,fig] = ...
    plotVirtualExperimentSegregationSequence(obj,axisDir,sequence,plotPoints,titleText)

  nRep = 200; %200

  % If no sequence is specified load one
  if(~exist('sequence') | isempty(sequence))
    sequence = obj.loadIterSequence();
  end
  
  if(~exist('plotPoints'))
    plotPoints = false;
  end

  if(isempty(sequence))
    disp('No sequence loaded')
    k = [];
    iter = [];
    fig = [];
    return
  end
  
  % If the first RetinalMap object in the sequence has no synapses,
  % then throw it out.
  if(numel(sequence) > 0 ...
     & sum(sequence(1).numPresynapticConnections) == 0)
    sequence = sequence(2:end);
    disp('Skipping first object in sequence, it has 0 synapses')
  end
  
  
  % Calculate the segregation at different distances for different
  % points in time.
  
  d = {};
  sep = {};
  k = zeros(numel(sequence),2);
  
  for i = 1:numel(sequence)
    switch(axisDir)
      case {'NT','AP','nt','ap'}
        [d{i},sep{i},k(i,:)] = ...
            sequence(i).virtualInjectionSegregationExperiment('MLfix',nRep);
      case {'DV','ML','dv','ml'}
        [d{i},sep{i},k(i,:)] = ...
            sequence(i).virtualInjectionSegregationExperiment('APfix',nRep);
    end
  end
  
  % Plot the figures
  
  if(obj.plotFigures == 0)
    disp('plotVirtualExperimentSegregationSequence: plotFigures = 0, hiding figures!')
    visFlag = 'off';
  else
    visFlag = 'on';
  end

  fig = figure('visible',visFlag);

  % Get the colour map
  col = winter(numel(sequence));
  col = col(end:-1:1,:); % We want it in opposite direction
  lineTypes = {'--','-'};
  
  hold on
  x = linspace(0,1,100);
  p = zeros(numel(sep),1);
  iter = zeros(numel(sep),1);  
  pLeg = {};
  
  if(plotPoints)
    for i = 1:numel(sep)
      plot(d{i},sep{i},'.','color',col(i,:));
    end
  end
  
  for i = 1:numel(sep)
    y = logistic(k(i,:),x);
    iter(i) = sequence(i).curStep/sequence(i).nSC;
    
    lt = lineTypes{mod(i,numel(lineTypes))+1};
    p(i) = plot(x,y,'color',col(i,:),'linestyle',lt,'linewidth',2);
    pLeg{i} = sprintf('Iter = %d', iter(i));
  end

  hold off
  legend(p,pLeg);
    
  ylabel('Segregation in Retina','fontsize',20) 

  switch(axisDir)
    case {'NT','AP','nt','ap'}
      xlabel('Anterior - Posterior separation','fontsize',20)
    case {'DV','ML','dv','ml'}
      xlabel('Medial - Lateral separation','fontsize',20)
  end

  set(gca,'fontsize',16)
  
  box off
  axis([0 1 0.5 1])
  
  if(exist('titleText'))
    title(titleText)
  end  
  
  fName = sprintf('%s/%s-virtual-exp-segregation-%s-sequence.pdf', ...
                  obj.figurePath, obj.simName, axisDir);
  
  switch(obj.plotFigures)
    case {0,1}
      saveas(gcf,fName,'pdf');
    case 2
      obj.addParametersToFigure(strrep(fName,'.pdf','.eps'));
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function y = logistic(k,x)
    % y = -0.5./(1+(x/inflection).^slope)+1;
    y = 1 - 0.5./(1+(x/k(1)).^k(2));
    
  end
    
end