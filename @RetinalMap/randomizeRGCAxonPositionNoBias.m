% !!!! Figure out how we do when we want initial bias, how do we project them
% One possible finding is that when we add initial bias RGC death decreases

function randomizeRGCAxonPositionNoBias(obj)

  disp('Placing initial axons without initial bias')

  mlRange = [min(obj.SCml) max(obj.SCml)];
  % apRange = [min(obj.SCap) max(obj.SCap)];

  obj.RGCml = mlRange(1) + (mlRange(2)-mlRange(1))*rand(obj.nRGC,1);

end