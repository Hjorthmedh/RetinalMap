function placeNeurons(obj)

  % I have removed the minNum parameter, because if we just naively
  % cull a dMin:ed population of positions, then the resulting
  % population will not have the right distribution of distances
  
  % Place Retinal neurons

  switch(obj.eyeType)
    case 'sphere'
      obj.placeRetinalRGCsphere();

      assert(strcmp(obj.gradientGenerationMethod,'diffusion'));

    case 'disk'
  
      obj.placeRetinalRGCdisk();
    
    otherwise
      fprintf('Unknown eye type %s (use sphere or disk)\n', obj.eyeType)
   
  end

  switch(obj.RGCdensity)
    case 'uniform'
      % Do nothing, we already have uniform density
      
    case 'varying'
      % Cull the RGC population to get relative densities matching Jerom 1998
      disp('!!! Please modify code to use placeRGCdiskRealDensity instead')
      beep
      keyboard
      % !!! This no longer works, since I am not using minNum
      % for naive retinal placement above
      obj.setRGCdensity();
  
    otherwise
      fprintf('placeNeurons: Unknown RGCdensity: %s\n', obj.RGCdensity)
      keyboard
  end
  
  % Place SC neurons
  obj.placeSC(obj.SC3Dflag);
  
  
end