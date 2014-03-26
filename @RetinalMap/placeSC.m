function placeSC(obj,make3D)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % 3D SC
    
  if(exist('make3D') & make3D == true)
    [obj.SCap,obj.SCml,obj.SCz, ...
     obj.SCapPadding,obj.SCmlPadding,obj.SCzPadding] = ...
        obj.dMin3D(@obj.filterSCPosition3D,obj.nSC);

    maxDim = max(max(obj.SCap)-min(obj.SCap),max(obj.SCml)-min(obj.SCml));
    minAP = min(obj.SCap);
    minML = min(obj.SCml);
    
    obj.SCap = (obj.SCap-minAP)/maxDim;
    obj.SCml = (obj.SCml-minML)/maxDim;
    obj.SCz = obj.SCz/maxDim;
    
    obj.SCapPadding = (obj.SCapPadding-minAP)/maxDim;
    obj.SCmlPadding = (obj.SCmlPadding-minML)/maxDim;
    obj.SCzPadding = obj.SCzPadding/maxDim;
    
    return
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % 2D SC
  
  % New default method, dMin
  [AP,ML,APpad,MLpad] = obj.dMin(@obj.filterSCPosition,obj.nSC);
    
  obj.SCap = AP;
  obj.SCml = ML;
  obj.SCapPadding = APpad;
  obj.SCmlPadding = MLpad;       
   
  % Rescale so that the largest dimension goes from 0 to 1
  maxDim = max(max(obj.SCap)-min(obj.SCap),max(obj.SCml)-min(obj.SCml));
  minAP = min(obj.SCap);
  minML = min(obj.SCml);
  
  obj.SCap = (obj.SCap-minAP)/maxDim;
  obj.SCml = (obj.SCml-minML)/maxDim;
  
  obj.SCapPadding = (obj.SCapPadding-minAP)/maxDim;
  obj.SCmlPadding = (obj.SCmlPadding-minML)/maxDim;
  
  [nnRI,vRI] = obj.nearestNeighbourStatistics('SC');
  
  fprintf('SC regularity indexes: nearest neighbour = %.3f, voronoi = %.3f\n', ...
          nnRI,vRI);
    
end
