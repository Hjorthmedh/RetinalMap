% This function takes Dans experimental fits for the different
% ages, and tries to find which 

function matchIterationToAgeSegregationExperiment(obj,injectionType,nRep)

  switch(injectionType)
    case {'MLfix','AP'}
      axisDir = 'AP';
    case {'APfix','ML'}
      axisDir = 'ML';
    otherwise
      fprintf('Unknown injectionType = %s\n', injectionType)
      keyboard
  end  
  
  [WT,B2KO,SCsize] = obj.getDanSegregationFits();

  switch(obj.phenotype)
    case 'WT'
      expData = WT;
    case 'Beta2KO'
      expData = B2KO;
    otherwise
      
      disp('Only have data for WT and B2KO')
      return
  end
  
  [age,bestIter] = findMatch(expData,SCsize);

  plotBestMatch(expData,SCsize,minIdx,kFit,bestIter/obj.nSC);

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function [age,bestIter] = findMatch(expData,SCsize)
  
    seq = obj.loadIterSequence();
    
    if(sum(seq(1).totalWeightRGC) == 0)
      % Discard t=0 if it has no synapses
      seq = seq(2:end);
    end
    
    fitError = zeros(numel(seq),numel(expData));
    kFit = zeros(numel(seq),2);
    
    for i = 1:numel(seq)
      [distance,segregation,kFit(i,:)] = ...
          seq(i).virtualInjectionSegregationExperiment(injectionType,nRep);
      
      for j = 1:numel(expData.age)
        % Calculate the error for each age
        s = logistic(expData.kAP(j,:), ...
                     scaleDistance(distance,SCsize,expData.age(j)));
        fitError(i,j) = sum((s-segregation).^2);
      end
    end
      
    [~,minIdx] = min(fitError);
    
    age = expData.age;
    bestIter = cat(1,seq(minIdx).curStep);
    
  end

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function plotBestMatch(expData,SCsize,minIdx,kFit,bestIter)
    
    if(obj.plotFigures == 0)
      disp('plotVirtualExperimentSegregation: plotFigures = 0, hiding figures!')
      visFlag = 'off';
    else
      visFlag = 'on';
    end

    fig = figure('visible',visFlag);
    
    % Palette from colorbrewer.com
    refColour = [255, 247, 236;
                 254, 232, 200;
                 253, 212, 158;
                 253, 187, 132;
                 252, 141, 89;
                 239, 101, 72;
                 215, 48, 31;
                 179, 0, 0;
                 127, 0, 0]/255;

    x = linspace(0,1,100);    
    for i = 1:numel(expData.age)
      yExp = logistic(expData.kAP(i,:), ...
                      scaleDistance(x,SCsize,expData.age(i)));
      yMod = logistic(kFit(minIdx(i),:),x);
      
      plot(x,yExp,'--','color',refColour(i,:))
      hold on
      p(i) = plot(x,yMod,'-','color',refColour(i,:),'linewidth',2);
      pLeg{i} = sprintf('P%d (%d)',expData.age(i),bestIter(i));
    end
    
    switch(axisDir)
      case {'NT','AP','nt','ap'}
        xlabel('Anterior - Posterior separation','fontsize',24)
      case {'DV','ML','dv','ml'}
        xlabel('Medial - Lateral separation','fontsize',24)
    end
    ylabel('Segregation in Retina','fontsize',24) 

    set(gca,'fontsize',20)
    
    box off
      
    
    legend(p,pLeg);
    
    fName = sprintf('%s/%s-virtual-exp-segregation-%s-match-iter-to-age.pdf', ...
                    obj.figurePath, obj.simName, axisDir);
    
    switch(obj.plotFigures)
      case {0,1}
        saveas(gcf,fName,'pdf');
      case 2
        obj.addParametersToFigure(strrep(fName,'.pdf','.eps'));
    end
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Takes a normalised distance, and scales it up to full size
  function dScaled = scaleDistance(d,SCsize,age)

    idx = find(SCsize(:,1) == age);
    switch(injectionType)
      case {'MLfix','AP'}
        dScaled = d*SCsize(idx,2);
      case {'APfix','ML'}
        dScaled = d*SCsize(idx,3);
      otherwise
        fprintf('Unknown injectionType = %s\n', injectionType)
        keyboard
    end
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

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end