% If d is not specified, it finds the closests SC neurons to the ML
% and AP coordinates specified.
%
% If d is specified, it finds all the SC neurons within d of the coordinate

function SCidx = findSCneuron(obj,ml,ap,d)

  if(~exist('d') | isempty(d))

    % Find the neuron(s) closest to the ML and AP coordinates
    SCidx = zeros(size(ml));

    for i = 1:numel(ml)
    
      [~,idx] = min(sqrt((obj.SCap-ap).^2+(obj.SCml-ml).^2));
      SCidx(i) = idx(1);
    
    end
    
  else
    
    % Do a virtual injection in the SC, find the SC neurons within
    % injection site at (ML,AP).
    
    assert(numel(ml) == 1 & numel(ap) == 1);
    
    SCidx = find((obj.SCap-ap).^2 + (obj.SCml-ml).^2 < d^2);
    
  end
  
end