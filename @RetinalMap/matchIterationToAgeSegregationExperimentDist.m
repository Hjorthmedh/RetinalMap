% This function fits a curve to virtual segregation experiments for
% each iteration, and sees which one fits best to the given data.

function [age,bestIter] =  matchIterationToAgeSegregationExperimentDist(obj,normaliseSize)

  if(~exist('normaliseSize'))
    normaliseSize = true;
  end
  
  % Internally we only work with normalised distances in the model.
  % If normaliseSize is used, then the experimental data is
  % normalised before we do the fits. However, if it is set to
  % false, the model fits are rescaled to micrometer, before the
  % fitness of the curves are calculated to the unscaled
  % experimental data.
  
  % If normaliseSize is not used, we assume that the SC in the
  % model is 2000 micrometer.
  SCmodelScale = 2000;
  
  % We use Dans parameters as starting points for nlinfit
  [WTfit,B2KOfit,SCsize] = obj.getDanSegregationFits();  
  
  % These files are for the AP axis
  %WTfile = 'segregation/segregation-WT-C57BL6J.csv';
  %B2KOfile = 'segregation/segregation-B2KO.csv';

  % Latest AP files from Dan Lyngholm, 2012-12-12
  WTfile = 'segregation/NN-C57BL6J-AP.csv';
  B2KOfile = 'segregation/NN-B2KO-AP.csv';
  
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
      [age,bestIter,bestIdx] = findMatchDistBATCH(WT);
    case 'B2KO'
      [age,bestIter,bestIdx] = findMatchDistBATCH(B2KO);
    otherwise
      fprintf('No data for phenotype %s\n', obj.phenotype)
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  % This function uses the same virtual fits, to compare with all age groups
  
  function [age,bestIter,bestIdx] = findMatchDistBATCH(expData)

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
    
    nInj = 100;
    
    fitScore = zeros(numel(seq),numel(expData));
    kFit = zeros(2,numel(seq));

    % Note to self: Could use squeeze to collapse dimensions, but choose to
    % reorder them instead.
    
    % Calculate all the fitting curves, nSets, each with nInj injections
    for i = 1:numel(seq)
      fprintf('Seq: %d/%d\n', i, numel(seq))
      
      % This is all done in normalised distance, always
      [distance,segregation,kFit(:,i)] = ...
          seq(i).virtualInjectionSegregationExperiment(injectionType,nInj);
    end
    
    disp('Calculating scores...')
    
    % For each timepoint of experiment, find the iteration that
    % matches best

    distTotal = zeros(numel(age),numel(seq));    
    
    for k = 1:numel(age)
      dataIdx = find(expData.age == age(k));
      
      for i = 1:numel(seq)
        score = zeros(numel(dataIdx),1);
        
        if(normaliseSize)
          % How far are the experimental data points from the virtual fit
          y = logistic(kFit(:,i),expData.distScaled(dataIdx));
          distTotal(k,i) = sum((expData.NNmean(dataIdx)-y).^2);
        else
          % Non-default behaviour, we want to compapre the unscaled
          % distances, however our fits are in normalised
          % distances, so we downscale distance by 2000 micrometers
          % regardless of age of the SC.
          
          y = logistic(kFit(:,i),expData.dist(dataIdx)/SCmodelScale);
          distTotal(k,i) = sum((expData.NNmean(dataIdx)-y).^2);
        end
      end
        
    end
    
    % Now find out for each age which iteration has smallest score
  
    xRange = linspace(0,1,100);
    
    for k = 1:numel(age)

      %Instead of naively using the score total, maybe do it with a
      %smoothed version of the score total instead... To account
      %for fluctuations.

      smoothWindow = 10; % With a reportStep of 20, this corresponds
                         % to 500 iteration smoothing window
      
      distTotalSmoothed = smooth(distTotal(k,:),smoothWindow);
      
      [~,idx] = min(distTotalSmoothed);

      bestIdx(k) = idx(1);      
      bestIter(k) = seq(idx(1)).curStep / seq(idx(1)).nSC;      
      
      if(1)
        % Debug plots
        
        figure
        subplot(2,1,1)
        
        hold on
        
        % Draw the fit line
        y = logistic(kFit(:,bestIdx(k)),xRange);
        dataIdx = find(expData.age == age(k));
        
        if(normaliseSize)
          plot(xRange,y,'k-');
          % Plot the data points
          plot(expData.distScaled(dataIdx),expData.NNmean(dataIdx), ...
               'k.', 'markersize', 20);

          xlabel('Normalised distance on SC','fontsize',20)
        else
          plot(xRange*SCmodelScale,y,'k-');
          
          % Plot the data points
          plot(expData.dist(dataIdx),expData.NNmean(dataIdx), ...
               'k.', 'markersize', 20);

          xlabel('SC distance','fontsize',20)
        end
        
        ylabel('Segregation in retina','fontsize',20)
        box off
        set(gca,'fontsize',15);
        title(sprintf('%s P%d (iter %d)', ...
                      expData.type, age(k), bestIter(k)))
        
        subplot(2,1,2)
        iter = cat(1,seq.curStep) ./ cat(1,seq.nSC);
        plot(iter,distTotal(k,:),'k-', ...
             iter,distTotalSmoothed,'r-', ...
             iter(bestIdx(k)),distTotalSmoothed(bestIdx(k)),'r*');
        xlabel('Iteration','fontsize',20)
        ylabel('Mismatch','fontsize',20)
        box off
        set(gca,'fontsize',15)
        
        fName = sprintf('%s/%s-map-age-to-iter-P%d-distMatch.pdf', ...
                  obj.figurePath, obj.simName, age(k));        
        saveas(gcf,fName,'pdf');

      end
      
    end
   catch e
    getReport(e)
    keyboard
   end

   try
     fName = sprintf('%s/%s-matchIterToAge-dist-save.mat', ...
                     obj.dataPath, obj.simName);
     save(fName,'distTotal','kFit','age','bestIdx','bestIter','smoothWindow','obj');
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