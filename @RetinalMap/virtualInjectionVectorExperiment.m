function [SCdistance,flipFlag] = virtualInjectionVectorExperiment(obj,axisDir,nRep)

  % !!! Consider changing the injections to be more systematic, so
  % that we do them at a range of fixed distances from each other, discard if
  % either are outside the retina.
  if(~exist('nRep'))
    nRep = 10000;
  end
  
  fprintf('Doing virtual injections for vector flip experiment\n')

  % Code currently only handles eye-disks, but easy to extend
  assert(strcmp(obj.eyeType,'disk'));

  % Plot the injections
  debugFlag = false;

  % Find all SC neurons within injection site, diameter approx 3.5%
  % of SC. Radius thus approximately 1.75% of SC.
  % !!! Update, from Dans thesis, chap 3.3.1 -- 2.8% of SC, diameter
  
  SCwidth = max(obj.SCap(:))-min(obj.SCap(:));  
  SCdist = SCwidth*0.028/2;
 
  flipFlag = zeros(nRep,1);
  SCdistance = zeros(nRep,1);
  
  for expCtr = 1:nRep
    % Select location of first SC injection
    okPoint = false;
  
    while(~okPoint)
      SC1ap = 0.3 + 0.4*rand(1);
      SC1ml = 0.3 + 0.4*rand(1);

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
      switch(axisDir)
        case {'ML','ml'} % AP fixed
          SC2ap = SC1ap;
          SC2ml = rand(1);
          
        case {'AP','ap'} % ML fixed
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

      if(~isempty(SCidx2))
        okPoint = 1;
      end
    end
    
    SCdistance(expCtr) = sqrt((SC1ap - SC2ap)^2 + (SC1ml - SC2ml)^2);

    % Find location of presynaptic RGC to the selected SC neurons
    [centNT1,centDV1] = calculateProjCentroid(SCidx1);
    [centNT2,centDV2] = calculateProjCentroid(SCidx2);    
    
    % Are they ordered correctly along the relevant axis?
    
    switch(axisDir)
      case {'AP','ap'}
        if((SC1ap < SC2ap & centNT1 > centNT2) ...
           | (SC1ap > SC2ap & centNT1 < centNT2))
          flipFlag(expCtr) = 0; % They are on the right side
        else
          flipFlag(expCtr) = 1;
        end
      
      case {'ML','ml'}
        if((SC1ml < SC2ml & centDV1 > centDV2) ...
           | (SC1ml > SC2ml & centDV1 < centDV2))
          flipFlag(expCtr) = 0;
        else
          flipFlag(expCtr) = 1;
        end
    end
    
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  function [centroidNT, centroidDV] = calculateProjCentroid(SCidx)

    if(isempty(SCidx))
      centroidNT = NaN;
      centroidDV = NaN;
      return
    end
    
    centroidNT = 0;
    centroidDV = 0;
    totalWeightSC = 0;
    
    for i = 1:numel(SCidx)
      idx = SCidx(i);
      totalWeightSC = totalWeightSC ...
          + sum(obj.presynapticWeight(1:obj.numPresynapticConnections(idx),idx));

      for j = 1:obj.numPresynapticConnections(idx)
        centroidNT = centroidNT ...
            + obj.presynapticWeight(j,idx) ...
              * obj.RGCnt(obj.presynapticConnections(j,idx));
        centroidDV = centroidDV ...
            + obj.presynapticWeight(j,idx) ...
              * obj.RGCdv(obj.presynapticConnections(j,idx));
      end
    end

    centroidNT = centroidNT/ totalWeightSC;
    centroidDV = centroidDV / totalWeightSC;
    
  end
  
  
end
