function setupGradients(obj)


  
  switch(obj.gradientGenerationMethod)
  
    case 'diffusion'
      
      % The diffusion calculation requires a certain minimum number
      % of points.
      assert(obj.nRGC >= 10000 & obj.nSC >= 10000);
      
      % Diffusion gradients in retina
      
      switch(obj.eyeType)
        case 'sphere'
          % 3D
          obj.setupGradientsRetina3Ddiffusion();
      
        case 'disk'
          % 2D
          obj.setupGradientsRetina2Ddiffusion();
      end
      
      % Diffusion gradients in SC
      obj.setupGradientsSCdiffusion();

      
    case 'naive'
      
      obj.makeNaiveGradients('withNoise');
      
    case 'phenotype'

      obj.loadGradients(obj.phenotype);
      
    otherwise

      fprintf('Non-default gradient generation method: %s\n', ...
              obj.gradientGenerationMethod)
      
      disp('Calling makeNaiveGradients.')
      
      obj.makeNaiveGradients(obj.gradientGenerationMethod,'withNoise');
      
  end
  
  % If kMask is nonzero, this function masks the gradients accordingly
  obj.maskGradients(1); 
  % 1 = plot resulting gradients, 0 = dont plot gradients
  % kMask determines type of masking, see comments above
  
  

end