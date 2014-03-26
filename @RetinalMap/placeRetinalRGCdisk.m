% 2D retina functions

%
% !!! update to use dMin function?
%

function placeRetinalRGCdisk(obj,minNum)
  % Create a disk centered around (0.5,0.5) of radie 0.5
  % filled with RGC neurons randomly placed
      
  % We need a relatively fine grid for accurate diffusion
  if(exist('minNum'))
    nRGC = max(obj.nRGC,minNum);
  else
    nRGC = obj.nRGC;
  end
    
  if(1)
    % Default method now uses dMin to place the neurons

    [NT,DV,NTpad,DVpad] = obj.dMin(@obj.filterRetinalPositionsDisk,obj.nRGC);
    
  else
    % Previous method used Halton sets
    disp('Using obsolete Halton set neuron placement')
    
    % Randomize out 2*nRGC points, remove those outside our valid region
    % then take the nRGC first points

    H = haltonset(2,'Skip',1e3,'Leap',1e2);
    H = scramble(H,'RR2');

    X0 = H(1:2*nRGC,:);
    % We multiply by 1.2 and subtract 0.1 to get from -0.1 to 1.1 
    % (ie a bit of padding around the disk)
    [NT,DV,NTpad,DVpad] = ...
        obj.filterRetinalPositionsDisk(1.2*X0(:,1)-0.1,1.2*X0(:,2)-0.1);

  end
    
  assert(length(NT) >= nRGC);

  obj.RGCnt = NT(1:nRGC);
  obj.RGCdv = DV(1:nRGC);
  obj.RGCntPadding = NTpad;
  obj.RGCdvPadding = DVpad;

  [nnRI,vRI] = obj.nearestNeighbourStatistics('retina');  
  
  fprintf('Retinal regularity index: nearest neighbour = %.3f, voronoi = %.3f\n', ...
          nnRI, vRI)
    
  if(~isempty(obj.nRGCtypes) & obj.nRGCtypes > 1)
    obj.RGCtype = ceil(obj.nRGCtypes * rand(obj.nRGC,1));
  end
  
end
