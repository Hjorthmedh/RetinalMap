% Set the RGC and SC colours according to a reference image

function [RGCcol,SCcol] = getNeuronColoursFromImage(obj,refImageName,flipImg)

  if(~exist('flipImg'))
    flipImg = 1;
  end

  % Load reference image, to colour the RGC neurons by
  img = imread(refImageName);
  
  % Make sure the image is square
  dim = min(size(img,1),size(img,2));
  img = img(1:dim,1:dim,:);
  
  RGCcol = zeros(obj.nRGC,3);
  SCcol = zeros(obj.nSC,3);
  
  X = linspace(max(obj.RGCnt),min(obj.RGCnt),size(img,2));
  if(flipImg)
    Y = linspace(max(obj.RGCdv),min(obj.RGCdv),size(img,1));    
  else
    Y = linspace(min(obj.RGCdv),max(obj.RGCdv),size(img,1));
  end
  
  for i = 1:3
    try
      RGCcol(:,i) = interp2(X,Y,double(img(:,:,i))/255.0,obj.RGCnt,obj.RGCdv);
    catch e
      getReport(e)
      keyboard
    end
    
  end
  
  % Find the correspoding colours for the SC neurons
  for i = 1:obj.nSC
      
    idx = obj.presynapticConnections(1:obj.numPresynapticConnections(i),i);
    w = obj.presynapticWeight(1:obj.numPresynapticConnections(i),i);
    
    for j = 1:numel(idx)
      SCcol(i,:) = SCcol(i,:) + RGCcol(idx(j),:)*w(j);
    end
    
    SCcol(i,:) = SCcol(i,:) / sum(w);
    
  end

 
  idx = find(isnan(SCcol(:,1)));
  SCcol(idx,:) = 0;

end