% This function simulates a virtual injection and calculates how
% large the resulting projection zone is
%

function [TZspreadInj,nRGCLabeled,nSCLabeled] = ...
      calculateTerminationZoneSizeInjection(obj,nRep,percentile,plotFlag)

  debugFigures = false; %true;
  nRGCLabeled = zeros(nRep,1);
  nSCLabeled = zeros(nRep,1);
  
  if(~exist('plotFlag'))
    plotFlag = false;
  end
  
  if(~exist('nRep') | isempty(nRep))
   nRep = 20; % Needed 100 for old method, lets try fewer for new.
  end 
  
  if(~exist('percentile'))
    percentile = 0.95;
  end
  
  % Find all SC neurons within injection site, diameter approx 3.5%
  % of SC. Radius thus approximately 1.75% of SC.
  % !!! Update, from Dans thesis, chap 3.3.1 -- 2.8% of SC, diameter
  % For WT adult, this should give 8.5% of retinal diameter
  
  SCwidth = max(obj.SCap(:))-min(obj.SCap(:));  
  SCdist = SCwidth*0.028/2;

  kEst = 0.03; % Kernel estimate, initial guess (negative sign in
               % front means use it as is)
  
  for iRep = 1:nRep
    
    % Randomize an injection site, make sure there is at least one
    % SC neuron close by
    injIdx = [];
    
    while(isempty(injIdx))
      injAP = 0.3+0.4*rand(1);
      injML = 0.3+0.4*rand(1);
      
      injIdx = find(sqrt((obj.SCap-injAP).^2 ...
                           + (obj.SCml-injML).^2) < SCdist);
    end
    
    fprintf('AP %.2f, ML %.2f: ', injAP, injML)
    [preIdx,weight] = obj.findPresynapticRGC(injIdx);
    
    nRGCLabeled(iRep) = numel(preIdx);
    nSCLabeled(iRep) = numel(injIdx);
    
    % The old way of doing it:
    %
    % [TZspreadInj(iRep),NTcentre,DVcentre] = projectionSpread(preIdx);
    %   
    % if(debugFigures & iRep == 1)
    %   plotInjectionProj(injIdx,preIdx,TZspreadInj(iRep),NTcentre,DVcentre)
    % end
    
    % The new way of doing things:
    
    % Note the kEst from previous iteration gets used as starting
    % point for next injection (to speed up convergence)
    approxFlag = true;
    %[TZspreadInj(iRep),kEst] = obj.contourAnalysis(preIdx,percentile,weight, ...
    %                                               injIdx,plotFlag,kEst,approxFlag);
    
    % Modified the code so it can handle multiple percentiles
    for i = 1:numel(percentile)
      if(i == 1)
        kEst = abs(kEst);
      else
        % Since we already estimated kEst for i=1, we can reuse it,
        % does not change with different percentiles, just with
        % different injections.
        kEst = -abs(kEst);
      end
        
      % The ones is because we assume all weights should be 1 for
      % the RGC.
      [TZspreadInj(iRep,i),kEst] = obj.contourAnalysisKDE(preIdx,percentile(i));
      
      %[TZspreadInj(iRep,i),kEst] = obj.contourAnalysis(preIdx,percentile(i),...
      %                                                 ones(size(preIdx)), injIdx, ...
      %                                                 plotFlag,kEst,approxFlag);
    end
    
  end

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Note that we are ignoring synaptic weights here, since the
  % experimental injections are only analysed binary, either there
  % is a connection or there is not.
  
  function [radius95,NTcentre,DVcentre] = projectionSpread(RGCidx)

    NTcentre = mean(obj.RGCnt(RGCidx));
    DVcentre = mean(obj.RGCdv(RGCidx));
    
    d = sqrt((obj.RGCnt(RGCidx)-NTcentre).^2 ...
             + (obj.RGCdv(RGCidx)-DVcentre).^2);
    
    % The radius that includes 95% of the projected RGC
    radius95 = prctile(d,95);
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function plotInjectionProj(SCidx,RGCidx,projRadius,NTproj,DVproj)

    try
    figure
    subplot(1,2,1)
    idxFreeRGC = setdiff(1:obj.nRGC,RGCidx);
    plot(obj.RGCnt(idxFreeRGC),obj.RGCdv(idxFreeRGC),...
         'o','color', [1 1 1]*0.8);
    hold on
    plot(obj.RGCnt(RGCidx),obj.RGCdv(RGCidx),'r*');
    fi = linspace(0,2*pi,100);
    x = projRadius*cos(fi)+NTproj;
    y = projRadius*sin(fi)+DVproj;
    plot(x,y,'k-');
    hold off
    xlabel('Temporal - Nasal')
    ylabel('Ventral - Dorsal')
    set(gca, 'xdir','reverse','ydir','reverse');
    axis equal
    axis([0 max(max(obj.RGCnt),1) 0 max(max(obj.RGCdv),1)])
    box off
    
    subplot(1,2,2)
    idxFreeSC = setdiff(1:obj.nSC,SCidx);
    plot(obj.SCap(idxFreeSC),obj.SCml(idxFreeSC), ...
         'o', 'color', [1 1 1]*0.8);
    hold on
    plot(obj.SCap(SCidx),obj.SCml(SCidx),'r*');
    hold off
    xlabel('Anterior - Posterior')
    ylabel('Medial - Lateral')
    axis equal
    axis([0 max(max(obj.SCap),1) 0 max(max(obj.SCml),1)])
    box off
    
    catch e
      getReport(e)
      keyboard
    end
    
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end