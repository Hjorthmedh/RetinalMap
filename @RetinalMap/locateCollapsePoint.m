function NT = locateCollapsePoint(obj,debugFlag)

if(~exist('debugFlag'))
  debugFlag = 0;
end

try

  % We want to locate the collapse point for the Isl2-EphA3 maps

  % For a small region of NT axis, find all synapses onto SC
  
  nBins = 50; %50;
  NTedges = linspace(min(obj.RGCnt),max(obj.RGCnt),nBins+1);
  twoMaps = zeros(nBins,1);
  
  nTot = sum(obj.numPresynapticConnections);
  synapseTargets = zeros(nTot,1);
  presynapticWeights = zeros(nTot,1);
  binMember = zeros(nTot,1);
  
  nextPos = 1;
  
  for i = 1:obj.nSC
    for j = 1:obj.numPresynapticConnections(i)
      idx = obj.presynapticConnections(j,i);
      w = obj.presynapticWeight(j,i);
 
      binIdx = find(NTedges(1:end-1) <= obj.RGCnt(idx) ...
                    & obj.RGCnt(idx) <= NTedges(2:end),1);
      
      binMember(nextPos) = binIdx;
      synapseTargets(nextPos) = i;
      presynapticWeights(nextPos) = w;
      
      nextPos = nextPos + 1;
      
    end
  end  
  
  % Perform K-means clustering on the AP location of the synapses
  %
  % This is also used in virtualInjectionMappingExperiment.m, if
  % you change here, update there accordingly.
  
  allMean = NaN*zeros(nBins,2);
  allStd = NaN*zeros(nBins,2);
  
  for i = 1:nBins

    binIdx = find(binMember == i);

    if(isempty(binIdx))
      twoMaps(i) = NaN;
    else
      SCap = obj.SCap(synapseTargets(binIdx));

      idx = kmeans(SCap,2,'emptyaction','singleton');
    
      idx1 = find(idx == 1);
      idx2 = find(idx == 2);
    
      n1 = numel(idx1);
      n2 = numel(idx2);
    
      nRatio = min(n1,n2)/max(n2,n2);
    
      mean1 = mean(SCap(idx1));
      std1  = std(SCap(idx1));
      mean2 = mean(SCap(idx2));
      std2  = std(SCap(idx2));
    
      allMean(i,1) = mean1;
      allMean(i,2) = mean2;
      allStd(i,1) = std1;
      allStd(i,2) = std2;
    
      separationReq = 1.5; % 2?
      minSepDist = 0.15;
      minRatio = 0.05; % This is seems to be pretty permissive
   
      % Unsure about criteria
      %if(abs(mean1-mean2) > minSepDist & nRatio > minRatio)
      if(abs(mean1-mean2) > (std1+std2)*separationReq & nRatio > minRatio) 
	twoMaps(i) = 1;
      else
	twoMaps(i) = 0;
      end
    end
  end
  
  % Find the point that separates the region with two projection
  % zones from that with one projection zone. Here the mean of the
  % two clusters must be separated by more than the mean of their
  % standard deviation for them to be considered two separate zones.
  
  if(debugFlag)

    obj.plotAxisMap('AP');
    
    [c1,idxM1] = max(allMean,[],2);
    [c2,idxM2] = min(allMean,[],2);
    s1 = allStd(sub2ind(size(allStd), ...
                        1:size(allStd,1),transpose(idxM1)));
    s2 = allStd(sub2ind(size(allStd), ...
                        1:size(allStd,1),transpose(idxM2)));
    
    nt = NTedges(1:end-1) + diff(NTedges)/2;
    figure, hold on
    errorbar(nt,c1,s1,'k');
    errorbar(nt,c2,s2,'r');
    
    tIdx = find(twoMaps);
    plot(nt(tIdx),c1(tIdx),'k.','markersize',20)
    plot(nt(tIdx),c2(tIdx),'r.','markersize',20)    
    
    hold off
    
  end
   
  % We exclude NaN values from the collapse point analysis

  if(nnz(isnan(twoMaps)) == numel(twoMaps))
    % Entire vector is NaN
    NT = NaN;
    return
  end

  %if(nnz(twoMaps) == 0)
  if(nnz(twoMaps) == nnz(isnan(twoMaps))) 
    % All are zero, or all non-zeros are NaN
    % Just one map
    NT = 0;
    return
  end
  
  if(nnz(twoMaps) == nBins)
    NT = inf;
    return
  end
  
  collapseIdx = find(twoMaps == 0,1,'first');
  
  NT = (NTedges(collapseIdx) + NTedges(collapseIdx+1))/2;

catch e
  getReport(e)
  keyboard
end
  
end
