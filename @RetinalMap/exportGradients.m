function exportGradients(obj)

  retinaFile = sprintf('%s/%s-retinal-gradients.txt', ...
                       obj.dataPath, obj.simName);
  SCFile = sprintf('%s/%s-SC-gradients.txt', ...
                   obj.dataPath, obj.simName);

  retMat = [obj.RGCnt, obj.RGCdv, ...
            obj.RGCEphA, obj.RGCEphB, ...
            obj.RGCephrinA, obj.RGCephrinB];
  
  SCmat = [obj.SCap, obj.SCml, ...
           obj.SCephrinA, obj.SCephrinB, ...
           obj.SCEphA, obj.SCEphB];
  
  fprintf('Saving retinal gradients to %s\n', retinaFile)
  disp('Format: nt dv EphA EphB ephrinA ephrinB')
  save(retinaFile,'retMat','-ascii','-double');
  fprintf('Saving SC gradients to %s\n', SCFile)
  disp('Format: ap ml ephrinA ephrinB EphA EphB')
  save(SCFile,'SCmat','-ascii','-double');

  obj.saveState();
  
end