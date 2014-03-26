function setNumberOfSC(obj,nSC)

  fprintf('Reducing the number of SC neurons from %d to %d\n', ...
          obj.nSC, nSC)

  assert(obj.nSC >= nSC);
  
  p = randperm(obj.nSC);
  idx = p(1:nSC);

  obj.nSC = nSC;

  obj.SCap = obj.SCap(idx);
  obj.SCml = obj.SCml(idx);
  obj.SCephrinA = obj.SCephrinA(idx);
  obj.SCephrinB = obj.SCephrinB(idx);

  if(~isempty(obj.SCEphA))
    obj.SCEphA = obj.SCEphA(idx);
  end

  if(~isempty(obj.SCEphB))
    obj.SCEphB = obj.SCEphB(idx);
  end


  % Clean the padding
  obj.SCapPadding = [];
  obj.SCmlPadding = [];
    
end