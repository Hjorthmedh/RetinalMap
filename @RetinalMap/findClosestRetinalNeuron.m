function closestIdx = findClosestRetinalNeuron(obj, X, Y, subsetIdx)

  if(~exist('subsetIdx'))
    subsetIdx = 1:obj.nRGC;
  end
  
  % Match dimensions
  
  if(numel(Y) == 1 & numel(X) > 1)
    Y = Y*ones(size(X));
  end
  
  if(numel(X) == 1 & numel(Y) > 1)
    X = X * ones(size(Y));
  end
  
  closestIdx = NaN*ones(size(X));    

  
  switch(obj.eyeType)
    case 'sphere'
      % We are in spherical coordinates
      theta = X;
      phi = Y;
    
      for i = 1:numel(theta)
        [~,idx] = min((obj.RGCtheta(subsetIdx) - theta(i)).^2 ...
                      + (obj.RGCphi(subsetIdx) - phi(i)).^2);

        % We only want to return one point
        closestIdx(i) = subsetIdx(idx(1));      
      end
  
    case 'disk'

      % The eye is a disk
      RGCnt = X;
      RGCdv = Y;
  
      for i = 1:numel(RGCnt)
        [~,idx] = min((obj.RGCnt(subsetIdx) - RGCnt(i)).^2 + (obj.RGCdv(subsetIdx) - RGCdv(i)).^2);

        % We only want to return one point
        closestIdx(i) = subsetIdx(idx(1));
      end
      
    otherwise

      fprintf('findClosestRetinalNeuron: Unknown eye type %s\n', obj.eyeType)
      keyboard
      
  end


 
  
end
