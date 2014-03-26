function obj = stepMarkerInduction(obj, nSteps) 
  dt = obj.info.dt;
  try
    for i=1:nSteps
      % Compute induced marker
      totalWeightSC = sum(obj.connectionMatrix, 2);
      obj.info.IA = (obj.connectionMatrix * obj.RGCEphA) ./ totalWeightSC;
      obj.info.IB = (obj.connectionMatrix * obj.RGCEphB) ./ totalWeightSC;
      % obj.totalWeightSC = totalWeightSC;

      % Update target marker
      obj.SCephrinA = obj.SCephrinA ...
          + dt*(obj.info.alpha*(1 - obj.info.RGCEphAScale.*obj.info.IA.*obj.SCephrinA) ...
                + obj.info.beta*laplacian('A'));

      obj.SCephrinB = obj.SCephrinB ...
          + dt*(obj.info.alpha * (obj.info.IB - obj.SCephrinB) ...
                + obj.info.beta*laplacian('B'));

      % Update synaptic weights
      Psi = ((obj.info.RGCEphAScale*obj.SCephrinA * transpose(obj.RGCEphA)) - 1).^2 ...
            + (ones(obj.nSC, 1)*transpose(obj.RGCEphB) ...
               - obj.SCephrinB*ones(1, obj.nRGC)).^2;

      Phi = exp(-Psi/(2*obj.info.kappa^2));

      S = obj.connectionMatrix + dt * obj.info.gamma * Phi;
      dS = S./(ones(obj.nSC, 1) * sum(S, 1)) - obj.connectionMatrix;
      obj.connectionMatrix = obj.connectionMatrix + dS;
      % disp(max(max(abs(dS))))
      % obj.connectionMatrix = S./(ones(obj.nSC, 1) * sum(S, 1));
      % obj.totalWeightRGC = transpose(sum(obj.connectionMatrix, 1));
      assert(~any(isnan(obj.SCephrinA)))
      obj.curStep = obj.curStep + 1;
    end

  catch e
    getReport(e)
    keyboard
  end

% !!! Note potential problem with symmetry of Laplacian: 
% !!! Conservation of mass.

function out = laplacian(type) 
  try
    switch(type) 
      case 'A'
        conc = obj.SCephrinA;
      case 'B'
        conc = obj.SCephrinB;
      otherwise 
        assert(false);
    end
    
    out = zeros(size(conc));
    for i=1:obj.nSC 
      out(i) = sum((conc(obj.neighbourSC{i}) - conc(i)));
      % ...
      % ./obj.info.neighbourSCdist{i});
      
    end
    
  catch e
    getReport(e)
    keyboard
  end
  
end 

end
