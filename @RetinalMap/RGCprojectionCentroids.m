function [RGCcentAP,RGCcentML] = RGCprojectionCentroids(obj)

  % Make sure the variable is initialised
  assert(~isempty(obj.presynapticConnections));
  
  nCon = zeros(obj.nRGC,1);
  RGCcentAP = zeros(obj.nRGC,1);
  RGCcentML = zeros(obj.nRGC,1);
  
  % Find the SC projection centroids of all the RGC
  for i = 1:obj.nSC
    for j = 1:obj.numPresynapticConnections(i)
      idx = obj.presynapticConnections(j,i);
      w = obj.presynapticWeight(j,i);
      RGCcentAP(idx) = RGCcentAP(idx) + w*obj.SCap(i);
      RGCcentML(idx) = RGCcentML(idx) + w*obj.SCml(i);
      nCon(idx) = nCon(idx) + w;
    end
  end
  
  % Internal check, these should match
  try
    assert(nnz(nCon-obj.totalWeightRGC) == 0);
  catch e
    getReport(e)
    keyboard
  end      
  
  RGCcentAP = RGCcentAP ./ nCon;
  RGCcentML = RGCcentML ./ nCon;
  
end