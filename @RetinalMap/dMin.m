% constrainFcn = @r.filterRetinalPositionsDisk, 
% @r.filterRetinalPositionSphere or 
%
% If dMin is not given a value will be estimated for it

function [X,Y,Xout,Yout] = dMin(obj,constrainFcn,nPoints,dMin,packingFactor)

  % The user should not specify both packing factor and dMin at the
  % same time
  if(exist('packingFactor') & ~isempty(packingFactor) ...
         & exist('dMin') & ~isempty(dMin))
    disp('dMin: Do not specify both dMin and packing factor at the same time')
    keyboard
  end
  
  nPlaced = 1;
  nRejected = 0;
  nMaxRejected = 1e3*nPoints;
  nOutside = 0;

  % If the user has specified a dMin value, then we do not try to
  % force the exact number of neurons placed in the area
  forceN = false;
  
  X = zeros(nPoints,1);
  Y = zeros(nPoints,1);

  Xout = zeros(nPoints,1);
  Yout = zeros(nPoints,1);
  
  % First seed point, to avoid have to handle empty case in the
  % computer intensive inner loop
  xy = rand(1000,2);
  [x,y] = constrainFcn(xy(:,1),xy(:,2));

  if(~exist('dMin') | isempty(dMin) | isnan(dMin))
    forceN = true;
    
    % Estimate area, to calculate suitable dMin
    [~,area] = convhull(x,y);  

    % optimalPackingFactor = 1/6*pi*sqrt(3) = 0.9069 -- Steinhaus 1999, p 202
    % best for dMin = 0.5473 -- Masharu Tanemura, 1979
    if(~exist('packingFactor') | isempty(packingFactor))
      packingFactor = obj.dMinPackingFactor;
    end
    
    r = sqrt(packingFactor*area/(nPoints*pi));
    dMin = 2*r;
    
    fprintf('Estimated area %.2f\n', area)
    fprintf('No dMin specified. Estimated dMin to %.4f\n', dMin)
  end
  
  X(1) = x(1);
  Y(1) = y(1);

  dMin2 = dMin^2;
  
  while(nPlaced < nPoints & nRejected < nMaxRejected)
    if(nPlaced > 1)
      fprintf('%d neurons placed (%d rejected)\n', nPlaced, nRejected)
    end
    
    % First get points that are within the allowed region
    xy = 1.2*rand(nPoints*10,2)-0.1;
    [~,~,~,~,mask] = constrainFcn(xy(:,1),xy(:,2));

    % Check the dMin criteria for each of the points
    for i = 1:size(xy,1)
      d = [(xy(i,1)-X(1:nPlaced)).^2 + (xy(i,2)-Y(1:nPlaced)).^2;
           (xy(i,1)-Xout(1:nOutside)).^2 + (xy(i,2)-Yout(1:nOutside)).^2];
      
      if(nPlaced < nPoints & d >= dMin2)
        if(mask(i))
          % Inside
          nPlaced = nPlaced + 1;
          X(nPlaced) = xy(i,1);
          Y(nPlaced) = xy(i,2);
        else
          % Outside
          nOutside = nOutside + 1;
          Xout(nOutside) = xy(i,1);
          Yout(nOutside) = xy(i,2);
        end
      else
        nRejected = nRejected + 1;
      end
      
    end
    
  end
  
  fprintf('%d neurons placed (%d outside, %d rejected)\n', ...
          nPlaced, nOutside, nRejected)
    
  X = X(1:nPlaced);
  Y = Y(1:nPlaced);
  
  Xout = Xout(1:nOutside);
  Yout = Yout(1:nOutside);

  % Check that enough points were placed
  if(nPlaced < nPoints & forceN)
    fprintf('We were unable to place %d points, retrying.\n', nPoints)
    [X,Y,Xout,Yout] = obj.dMin(constrainFcn,nPoints);
  end
  
  % plot(X,Y,'k.')
  
end