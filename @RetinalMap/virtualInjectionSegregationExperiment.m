% This code requires that there are at least five points for each
% injection site, otherwise NaN is returned
%
% SC -> Retina
%
% injectionType: 'APfix' or 'MLfix'
% nRep: number of repetitions
% maxDist: largest allowed distance between injections (for later
% iterations we can ignore higher distances).
% k0: In case you want to seed the nlinfit startingpoint

function [distance, segregation, k, kLow, kHigh] = ...
  virtualInjectionSegregationExperiment(obj, injectionType, nRep, maxDist, k0)

  if(~exist('maxDist') | isempty(maxDist))
    maxDist = inf;
  end
  
  if(~exist('k0'))
    k0 = [];
  end
  
  fprintf('Doing virtual injection (%s), %d experiments\n', ...
          injectionType, nRep)
  
  % Code currently only handles eye-disks, but easy to extend
  assert(strcmp(obj.eyeType,'disk'));

  % Plot the injections
  debugFlag = false;

  % Find all SC neurons within injection site, diameter approx 3.5%
  % of SC. Radius thus approximately 1.75% of SC.
  % !!! Update, from Dans thesis, chap 3.3.1 -- 2.8% of SC, diameter
  
  SCwidth = max(obj.SCap(:))-min(obj.SCap(:));  
  SCdist = SCwidth*0.028/2;
  
  % Pre allocate
  segregation = zeros(nRep,1);
  distance = zeros(nRep,1);
  
  for expCtr = 1:nRep
  
    % [SC1ap, SC1ml, SC2ap, SC2ml] = getTwoPointsNaive(injectionType);
    [SC1ap, SC1ml, SC2ap, SC2ml] = getTwoPointsEvenSpacing(injectionType);    

    distance(expCtr) = sqrt((SC1ap - SC2ap)^2 + (SC1ml - SC2ml)^2);
  
    % Find location of presynaptic RGC to the selected SC neurons
    RGCidx1 = obj.findPresynapticRGC(SCidx1);
    RGCidx2 = obj.findPresynapticRGC(SCidx2);
  
    % Calculate segregation
    segregation(expCtr) = calculateSegregation(RGCidx1, RGCidx2);
    
    % if(isnan(segregation(expCtr)))
    %   disp('Got NaN segregation value')
    %   keyboard
    % end
    
    if(debugFlag)
      plotVirtualInjection(SCidx1,RGCidx1,SCidx2,RGCidx2);
    
      fprintf('SC: AP %.3f ML %.3f, SC: AP: %.3f ML: %.3f (d=%.3f,sep=%.3f)\n', ...
              SC1ap, SC1ml, SC2ap, SC2ml, distance(expCtr),segregation(expCtr))
    
      %keyboard
    end  
    
  end 
  
  % For Math5 there are some injections with SC neurons that lack
  % presynaptic RGC, remove those.
  
  idx = find(~isnan(segregation));
  if(numel(idx) < numel(segregation))
    fprintf('Warning, %d injected SC neurons had no presynatpic RGC\n', ...
            numel(segregation)-numel(idx))
  end
  
  segregation = segregation(idx);
  distance = distance(idx);
  
  if(numel(idx) < 10)
    fprintf('virtualInjectionSegregationExperiment: Only %d points valid\n', numel(idx))
    disp('No logistic curve fitted.')
    k = [NaN NaN];
  else
    [k,kLow,kHigh] = obj.fitSegregation(distance,segregation,k0);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  function [selfCloser, otherCloser] = ...
      countClosestNeighbours(RGCidx1, RGCidx2)

    selfCloser = 0;
    otherCloser = 0;
    
    % Handle Math5 case when there might not be enough neurons
    % --- calculateSegregation willl no longer allow empty cases, so
    % this should not occur anymore...
    if(isempty(RGCidx1))
      selfCloser = 0;
      otherCloser = 0;
      return
    elseif(isempty(RGCidx2))
      selfCloser = numel(RGCidx1);
      otherCloser = 0;
      return
    end

    for i = 1:numel(RGCidx1)
      
      if(ismember(RGCidx1(i),RGCidx2))
        % RGC is targeting both injection sites
        selfCloser = selfCloser + 0.5;
        otherCloser = otherCloser + 0.5;
      else
        % The neuron we are looking at is not double labeled
        
        % Compare to all in its own group, but itself
        nsIdx = setdiff(RGCidx1,RGCidx1(i));
      
        dOwn = min(sqrt((obj.RGCnt(RGCidx1(i)) - obj.RGCnt(nsIdx)).^2 ...
                        + (obj.RGCdv(RGCidx1(i)) - obj.RGCdv(nsIdx)).^2));
      
        % Not a member of RGCidx2, we already checked that
        dOther = min(sqrt((obj.RGCnt(RGCidx1(i)) - obj.RGCnt(RGCidx2)).^2 ...
                          + (obj.RGCdv(RGCidx1(i)) - obj.RGCdv(RGCidx2)).^2));
      
        try
          if(dOwn < dOther)
            selfCloser = selfCloser + 1;
          elseif(dOwn > dOther)
            otherCloser = otherCloser + 1;
          else
            selfCloser = selfCloser + 0.5;
            otherCloser = otherCloser + 0.5;
          end
        catch e
          getReport(e)
          keyboard
        end
      end
    end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function sep = calculateSegregation(RGCidxA, RGCidxB)
  
    if(isempty(RGCidxA) | isempty(RGCidxB))
    % if(numel(RGCidxA) < 5 | numel(RGCidxB) < 5)
      sep = NaN;
      return
    end
    
    [selfCloserA,otherCloserA] = countClosestNeighbours(RGCidxA, RGCidxB);
    [selfCloserB,otherCloserB] = countClosestNeighbours(RGCidxB, RGCidxA);
     
    sep = (selfCloserA + selfCloserB) ...
        / (selfCloserA + selfCloserB + otherCloserA + otherCloserB);
    
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
  
  function [SC1ap, SC1ml, SC2ap, SC2ml] = getTwoPointsNaive(injectionType)
  
    % Select location of first SC injection
    okPoint = false;
  
    while(~okPoint)
      
      switch(injectionType)
        case 'APfix'
          SC1ap = 0.3 + 0.4*rand(1); % Place it fairly centrally
          SC1ml = rand(1);
        
        case 'MLfix'
          SC1ap = rand(1);
          SC1ml = 0.3 + 0.4*rand(1); % Place it fairly centrally
          
        otherwise
          fprintf('Unknown injectionType: %s\n', injectionType)
          keyboard
      end

      % Make sure there is at least one SC neuron close to the
      % injection site
      SCidx1 = find(sqrt((obj.SCap - SC1ap).^2 + (obj.SCml - SC1ml).^2) ...
                    < SCdist);

      if(~isempty(SCidx1))
        okPoint = 1;
      end
      
    end
  
    % Second injection site, this one is along AP or ML axis
    okPoint = false;

    while(~okPoint)
      switch(injectionType)
        case 'APfix'
          SC2ap = SC1ap;
          SC2ml = rand(1);
          
        case 'MLfix'
          SC2ap = rand(1);
          SC2ml = SC1ml;
          
        otherwise
          fprintf('Unknown injectionType: %s\n', injectionType)
          keyboard
      end

      % Make sure there is at least one SC neuron close to the
      % injection site
      SCidx2 = find(sqrt((obj.SCap - SC2ap).^2 + (obj.SCml - SC2ml).^2) ...
                    < SCdist);

      injDist = sqrt((SC1ap-SC2ap)^2 + (SC1ml-SC2ml)^2);
      
      if(~isempty(SCidx2) & injDist < maxDist)
        okPoint = 1;
      end
    end
     
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % This function tries to make sure all injection distances are
  % equally represented
  
  function [SC1ap, SC1ml, SC2ap, SC2ml] = getTwoPointsEvenSpacing(injectionType)
  
    switch(injectionType)
      case 'APfix'
        dMax = max(obj.SCml) - min(obj.SCml); 
      case 'MLfix'
        dMax = max(obj.SCap) - min(obj.SCap);         
      otherwise
        fprintf('Unknown injectionType: %s\n', injectionType)
        keyboard
    end
    
    % If we used the max distance as max, then we would only have
    % one pair there.
    d = 0.9*dMax*rand(1);
    
    % The user can specify a maxdist also
    d = min(d,maxDist);
    
    nTrials = 0;
    maxTrials = 1000;
    
    okPoint = false;
    
    while(~okPoint)
      if(nTrials > maxTrials)
        disp(['virtualInjectionSegregationExperiment: '...
              'Maximum number of trials reached. Something is fishy...'])
        keyboard
      end
      
      % Find two points d apart
      
      switch(injectionType)
        case 'APfix'
          SC1ap = 0.3 + 0.4*rand(1);
          SC2ap = SC1ap;
          
          SC1ml = (1-d)*rand(1);
          SC2ml = SC1ml+d;
          
        case 'MLfix'
          SC1ml = 0.3 + 0.4*rand(1);
          SC2ml = SC1ml;
          
          SC1ap = (1-d)*rand(1);
          SC2ap = SC1ap+d;
      end      
      
      % Make sure there are SC neurons close to both points
      SCidx1 = find(sqrt((obj.SCap - SC1ap).^2 + (obj.SCml - SC1ml).^2) ...
                    < SCdist);      
      
      SCidx2 = find(sqrt((obj.SCap - SC2ap).^2 + (obj.SCml - SC2ml).^2) ...
                    < SCdist);

      if(~isempty(SCidx1) & ~isempty(SCidx2))
        okPoint = true;
      end
      
      nTrials = nTrials + 1;
    end
    
  end
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
end
