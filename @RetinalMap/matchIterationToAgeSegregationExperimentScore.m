function [age,bestIter] =  matchIterationToAgeSegregationExperimentScore(obj)

  disp('Deprecated, use matchIterationToAgeSegregationExperimentDist.m')
  assert(false)
  
  % We use Dans parameters as starting points for nlinfit
  [WTfit,B2KOfit,SCsize] = obj.getDanSegregationFits();  
  
  % These files are for the AP axis
  WTfile = 'segregation/segregation-WT-C57BL6J.csv';
  B2KOfile = 'segregation/segregation-B2KO.csv';
  
  injectionType = 'MLfix';
  
  WTdata = importdata(WTfile);
  B2KOdata = importdata(B2KOfile);
  
  WT.age = WTdata.data(:,1);
  WT.dist = WTdata.data(:,2);
  WT.distScaled = scaleSCdist(WT.dist,WT.age);
  WT.NNred = WTdata.data(:,3)/100;
  WT.NNgreen = WTdata.data(:,4)/100;
  WT.NNmean = WTdata.data(:,5)/100;
  WT.type = 'WT';

  B2KO.age = B2KOdata.data(:,1);
  B2KO.dist = B2KOdata.data(:,2);
  B2KO.distScaled = scaleSCdist(B2KO.dist,B2KO.age);
  B2KO.NNred = B2KOdata.data(:,3)/100;
  B2KO.NNgreen = B2KOdata.data(:,4)/100;
  B2KO.NNmean = B2KOdata.data(:,5)/100;
  B2KO.type = 'B2KO'
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  switch(obj.phenotype)
    case {'WT','wt'}
      [age,bestIter,bestIdx] = findMatchScore(WT);
    case 'B2KO'
      [age,bestIter,bestIdx] = findMatchScore(B2KO);
    otherwise
      fprintf('No data for phenotype %s\n', obj.phenotype)
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % This function uses different virtual experiments for the
  % iterations for different age groups, to better replicate the real experiments
  
  function [age,bestIter,bestIdx] = findMatchScore(expData)

    age = unique(expData.age);
    bestIter = zeros(size(age));
    bestIdx = zeros(size(age));

    skipTable = 2; % This clears more memory
    seq = obj.loadIterSequence(skipTable);

    if(sum(seq(1).totalWeightRGC) == 0)
      % Discard t=0 if it has no synapses
      seq = seq(2:end);
    end
    
    %matlabpool open
    
    for ia = 1:numel(age)
      % We want to seed the nlinfit with an age appropriate k0
      % value. Take Dans original fittings as start points
      switch(obj.phenotype)
        case {'WT','wt'}
          idxFitWT = find(WTfit.age == age(ia));
          k0 = WTfit.kAPscaled(idxFitWT,:);
        case {'B2KO'}
          idxFitB2KO = find(B2KOfit.age == age(ia));
          k0 = B2KOfit.kAPscaled(idxFitB2KO,:);
        otherwise
          fprintf('Iteration to age matching: Unsupported phenotype: %s\n',...
                  obj.phenotype)
          keyboard
      end

      [bestIter(ia),bestIdx(ia)] = findMatchScoreAge(expData,age(ia),seq,k0);
    end
    
    % matlabpool close
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
  % Helper function to findMatchScore above
  
  function [bestGroupIter,bestGroupIdx] = findMatchScoreAge(expData,groupAge,sequence,k0)

    nInj = 20; %10;
    nSets = 21; % 11;
    
    smoothWindow = 25; % With a reportStep of 20, this corresponds
                       % to 500 iteration smoothing window
    
    % Find maxDist from experimental data
    dataIdx = find(expData.age == groupAge);
    
    % The 0.2 comes from the fact that if we happen to be below the
    % inflection point with all points, we are going to make baad
    % fits (got that problem for WT P16)
    maxDist = min(1,max(expData.distScaled(dataIdx))+0.2);
      
    fprintf('For %s age P%d, using maxDist = %f.3\n',expData.type,groupAge,maxDist)
    
    
    fitScore = zeros(numel(sequence),1);
    kFit = zeros(nSets,2,numel(sequence));
       
    
    % Calculate all the fitting curves, nSets, each with nInj injections
    %parfor i = 1:numel(sequence)
    for i = 1:numel(sequence)
      fprintf('Seq: %d/%d\n', i, numel(sequence))
      for j = 1:nSets
        [distance,segregation,kFit(j,:,i)] = ...
            sequence(i).virtualInjectionSegregationExperiment(injectionType,nInj,maxDist);
      end
    end
    
    disp('Calculating scores...')

    % For each timepoint of experiment, find the iteration that
    % matches best

    scoreTotal = zeros(numel(sequence),1);    
    scoreTotalSmoothed = zeros(numel(sequence),1);
    
    for i = 1:numel(sequence)
      
      score = zeros(numel(dataIdx),1);
        
      for m = 1:numel(dataIdx)
          
        yAll = logisticALT(kFit(:,:,i),expData.distScaled(dataIdx(m)));
        
        % The score for that particular point is defined as the
        % number of curves you need to cross (counting from the
        % median curve), until you get to the data point.
          
        yAll = sort([inf; yAll; -inf]);
          
        idx = find(yAll(1:end-1) < expData.NNmean(dataIdx(m)) ...
                   & expData.NNmean(dataIdx(m)) <= yAll(2:end));
               
        assert(~isempty(idx));
          
        score(m) = abs(idx - ceil(nSets/2) - 1);
      end
        
      scoreTotal(i) = sum(score);
        
    end
    
    % Smooth the score total
    scoreTotalSmoothed = smooth(scoreTotal,smoothWindow);
    % scoreTotalSmoothed = smooth(scoreTotal,'loess',smoothWindow);    
       
    
    xRange = linspace(0,maxDist,100);

    % Find out which iteration has the lowest score
    
    %Instead of naively using the score total, maybe do it with a
    %smoothed version of the score total instead... To account
    %for fluctuations.
    % [~,idx] = min(scoreTotal);
    [~,idx] = min(scoreTotalSmoothed);    

    bestGroupIdx = idx(1);      
    bestGroupIter = sequence(idx(1)).curStep / sequence(idx(1)).nSC;

    try
      fName = sprintf('%s/%s-matchIterToAge-save-P%d.mat', ...
                      obj.dataPath, obj.simName,groupAge);
      save(fName,'scoreTotal','scoreTotalSmoothed','bestGroupIdx','bestGroupIter','obj');
    catch e
      getReport(e)
      keyboard
    end
    
    if(1)
      % Debug plots
        
      figure
      subplot(2,1,1)
        
      hold on
        
      % Draw the fit lines
      for j = 1:nSets
        y = logistic(kFit(j,:,bestGroupIdx),xRange);
        plot(xRange,y,'k-');
      end

      % !!! Only plot as far as we got data
      
      % Plot the data points
      plot(expData.distScaled(dataIdx),expData.NNmean(dataIdx), ...
           'k.', 'markersize', 20);
      xlabel('Distance in SC','fontsize',20)
      ylabel('Segregation in retina','fontsize',20)
      set(gca,'fontsize',15);
      title(sprintf('%s P%d (iter %d)', ...
                    expData.type, groupAge, bestGroupIter))
        
      subplot(2,1,2)
      iter = cat(1,sequence.curStep) ./ cat(1,sequence.nSC);
      plot(iter,scoreTotal/numel(dataIdx),'k-', ...
           iter,scoreTotalSmoothed/numel(dataIdx),'r-', ...
           iter(bestGroupIdx),scoreTotalSmoothed(bestGroupIdx)/numel(dataIdx),'r*');
      xlabel('Iteration','fontsize',20)
      ylabel('Average score','fontsize',20)
      set(gca,'fontsize',15)
      %axis([0 1 0 1])  
      
      fName = sprintf('%s/%s-map-age-to-iter-P%d.pdf', ...
                      obj.figurePath, obj.simName, groupAge);        
      saveas(gcf,fName,'pdf');
    
    end
    
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % This function uses the same virtual fits, to compare with all age groups
  
  function [age,bestIter,bestIdx] = findMatchScoreBATCH(expData)

  try
    age = unique(expData.age);
    bestIter = zeros(size(age));
    bestIdx = zeros(size(age));
    
    skipTable = 2; % Conserve more memory
    seq = obj.loadIterSequence(skipTable);
    
    if(sum(seq(1).totalWeightRGC) == 0)
      % Discard t=0 if it has no synapses
      seq = seq(2:end);
    end
    
    nInj = 30; %20
    nSets = 11; %101;
    
    fitScore = zeros(numel(seq),numel(expData));
    kFit = zeros(nSets,2,numel(seq));

    % Note to self: Could use squeeze to collapse dimensions, but choose to
    % reorder them instead.
    
    % Calculate all the fitting curves, nSets, each with nInj injections
    for i = 1:numel(seq)
      fprintf('Seq: %d/%d\n', i, numel(seq))
      for j = 1:nSets
        [distance,segregation,kFit(j,:,i)] = ...
            seq(i).virtualInjectionSegregationExperiment(injectionType,nInj);
      end
    end
    
    disp('Calculating scores...')
    
    % For each timepoint of experiment, find the iteration that
    % matches best

    scoreTotal = zeros(numel(age),numel(seq),1);    
    
    for k = 1:numel(age)
      dataIdx = find(expData.age == age(k));
      
      for i = 1:numel(seq)
        score = zeros(numel(dataIdx),1);
        
        for m = 1:numel(dataIdx)
          
          yAll = logisticALT(kFit(:,:,i),expData.distScaled(dataIdx(m)));
        
          % The score for that particular point is defined as the
          % number of curves you need to cross (counting from the
          % median curve), until you get to the data point.
          
          yAll = sort([inf; yAll; -inf]);
          
          idx = find(yAll(1:end-1) < expData.NNmean(dataIdx(m)) ...
                     & expData.NNmean(dataIdx(m)) <= yAll(2:end));
               
          assert(~isempty(idx));
          
          score(m) = abs(idx - ceil(nSets/2) - 1);
        end
        
        scoreTotal(k,i) = sum(score);
        
      end
        
    end
    
    % Now find out for each age which iteration has smallest score
  
    xRange = linspace(0,1,100);
    
    for k = 1:numel(age)

      %Instead of naively using the score total, maybe do it with a
      %smoothed version of the score total instead... To account
      %for fluctuations.
      [~,idx] = min(scoreTotal(k,:));

      bestIdx(k) = idx(1);      
      bestIter(k) = seq(idx(1)).curStep / seq(idx(1)).nSC;

      
      if(1)
        % Debug plots
        
        figure
        subplot(2,1,1)
        
        hold on
        
        % Draw the fit lines
        for j = 1:nSets
          y = logistic(kFit(j,:,bestIdx(k)),xRange);
          plot(xRange,y,'k-');
        end

        dataIdx = find(expData.age == age(k));
        
        % Plot the data points
        plot(expData.distScaled(dataIdx),expData.NNmean(dataIdx), ...
             'k.', 'markersize', 20);
        xlabel('Distance in SC','fontsize',20)
        ylabel('Segregation in retina','fontsize',20)
        set(gca,'fontsize',15);
        title(sprintf('%s P%d (iter %d)', ...
                      expData.type, age(k), bestIter(k)))
        
        subplot(2,1,2)
        iter = cat(1,seq.curStep) ./ cat(1,seq.nSC);
        plot(iter,scoreTotal(k,:),'k-', ...
             iter(bestIdx(k)),scoreTotal(k,bestIdx(k)),'r*');
        xlabel('Iteration','fontsize',20)
        ylabel('Total score','fontsize',20)
        set(gca,'fontsize',15)
        
        fName = sprintf('%s/%s-map-age-to-iter-P%d.pdf', ...
                  obj.figurePath, obj.simName, age(k));        
        saveas(gcf,fName,'pdf');

        
      end
      
    end
   catch e
    getReport(e)
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

  % This function assumes x and y are scalar
  
  function y = logisticALT(k,x)
    y = 1 - 0.5./(1+(x./k(:,1)).^k(:,2));
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Daniel Lyngholm, email: 18 july 2012
  % Age		AP		ML
  % 0		1910	1800
  % 2		1920	1960
  % 4		2000	1900
  % 6		1980	1850
  % 8		1890	1840
  % 12		1920	1890
  % 16		1970	1980
  % 22		1970	1960
  % 60		2072	2008  
  
  function dScaled = scaleSCdist(dRaw,age)
    % We want the largest axis to be between 0 and 1, so scale by AP size 
  
    SCinfo = [0, 1910, 1800;
              2, 1920, 1960;
              4, 2000, 1900;
              6, 1980, 1850;
              8, 1890, 1840;
              12, 1920, 1890;
              16, 1970, 1980;
              22, 1970, 1960;
              60, 2072, 2008];
        
    SCage = SCinfo(:,1);
    SCAP = SCinfo(:,2);
    SCML = SCinfo(:,3);
    
    scaleFactor = interp1(SCage,SCAP,age);
    
    dScaled = dRaw ./ scaleFactor;
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end