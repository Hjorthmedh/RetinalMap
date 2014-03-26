% Internal helper function
%
% This returns all the presynaptic RGC to the SC neurons specified.
% No information about which one is connected to which is returned
% here. This function is useful when doing virtual retrograde injections

function [preIdx,weight] = findPresynapticRGC(obj, SCidx)

  % % Pre-allocate to speed things up
  %
  % nMax = sum(obj.numPresynapticConnections(SCidx));
  % nextIdx = 1;
  % preIdx = zeros(nMax,1);
  %  
  % for i = 1:numel(SCidx)
  %   nCon = obj.numPresynapticConnections(SCidx(i));
  %  
  %   preIdx(nextIdx:(nextIdx+nCon-1)) = ...
  %       obj.presynapticConnections(1:nCon,SCidx(i));
  %  
  %   nextIdx = nextIdx + nCon;
  % end
  %
  % preIdx = unique(preIdx);
    
  % Calculate the weight for each RGC
  RGCweight = zeros(obj.nRGC,1);
  
  for i = 1:numel(SCidx)
    nCon = obj.numPresynapticConnections(SCidx(i));
    idx = obj.presynapticConnections(1:nCon,SCidx(i));
    RGCweight(idx) = RGCweight(idx) + obj.presynapticWeight(1:nCon,SCidx(i));
  end

  preIdx = find(RGCweight);
  weight = RGCweight(preIdx);
  
end
