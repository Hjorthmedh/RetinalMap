% idx is optional, if not given a random set of RGC are kept
function setNumberOfRGC(obj,nRGC,idx)

  fprintf('Reducing the number of RGC from %d to %d\n', ...
          numel(obj.RGCnt), nRGC)

  assert(numel(obj.RGCnt) >= nRGC);

  if(exist('idx'))
    % The user specified which RGC to keep
    % check integrity of inparameters
    assert(numel(idx) == nRGC);
  else
    % Keep only a random selection
    p = randperm(obj.nRGC);
    idx = p(1:nRGC);
  end

  switch(obj.eyeType)
    case 'sphere'
      obj.RGCtheta = obj.RGCtheta(idx);
      obj.RGCphi = obj.RGCphi(idx);
    case 'disk'
      obj.RGCnt = obj.RGCnt(idx);
      obj.RGCdv = obj.RGCdv(idx);
    otherwise
      fprintf('setNumbeOfRGC: Unknown eye type %s\n', obj.eyeType)
      keyboard
  end
  
  nRGCold = obj.nRGC;
  obj.nRGC = nRGC;
  
  obj.RGCEphA = obj.RGCEphA(idx);
  obj.RGCEphB = obj.RGCEphB(idx);
  
  if(~isempty(obj.RGCephrinA))
    obj.RGCephrinA = obj.RGCephrinA(idx);
  end
  
  if(~isempty(obj.RGCephrinA))
    obj.RGCephrinB = obj.RGCephrinB(idx);
  end

  if(~isempty(obj.Isl2PositiveRGC))
    % idx is a random permutation of a subset of the RGC
    mask = zeros(nRGCold,1);
    mask(obj.Isl2PositiveRGC) = 1;
    obj.Isl2PositiveRGC = find(mask(idx));
  end
  
  % Clean the padding
  obj.paddingPhi = [];
  obj.paddingTheta = [];

  obj.RGCntPadding = [];
  obj.RGCdvPadding = [];

end