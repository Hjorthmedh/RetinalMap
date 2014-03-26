function placeInitialSynapses(obj)

  if(~obj.useLocalJumps)
    disp('No initial synapses for global jump simulations')
    return
  end
      
  fprintf('Placing %d initial synapses per RGC axon\n', ...
          obj.nInitSynapses)
  
  rWidth = obj.RGCwidth * (max(obj.SCml) - min(obj.SCml));

  for RGCidx = 1:obj.nRGC

    SCidx = find(obj.RGCml(RGCidx) - rWidth/2 <= obj.SCml ...
                 & obj.SCml <= obj.RGCml(RGCidx) + rWidth/2);
        
    SCidx = SCidx(randperm(numel(SCidx)));
    SCidx = SCidx(1:min(obj.nInitSynapses,numel(SCidx)));

    for i = 1:numel(SCidx)
      obj.addSynapse(SCidx(i),RGCidx);
    end
  end
  
end