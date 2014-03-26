function makeMapMovie(obj,iterSequence, RGCidx, SCidx, resolution, frameRate)

  if(~exist('iterSequence') | isempty(iterSequence))
    % Load the iter sequence of the current object.
    iterSequence = obj.loadIterSequence();
  end
  
  if(~exist('RGCidx'))
    RGCidx = [];
  end
  
  if(~exist('SCidx'))
    SCidx = [];
  end
  
  if(~exist('resolution'))
    resolution = 500;
  end
  
  if(~exist('frameRate'))
    frameRate = 5;
  end
  
  assert(strcmp(iterSequence(1).eyeType,'disk'))
    
  if(~exist('MOVIE'))
    mkdir('MOVIE');
  end
  
  tic
  
  if(~isempty(SCidx))
    SCtags = sprintf('-%d',SCidx);
    SCtags = sprintf('-SC%s',SCtags);
  else
    SCtags = [];
  end
  
  if(~isempty(RGCidx))
    RGCtags = sprintf('-%d',RGCidx);
    RGCtags = sprintf('-RGC%s', RGCtags);
  else
    RGCtags = [];
  end
  
  movName = sprintf('MOVIE/%s-map-development%s%s.avi', ...
                    iterSequence(1).simName,RGCtags,SCtags);  
   
  vidObj = VideoWriter(movName);
  vidObj.Quality = 100;
  vidObj.FrameRate = frameRate;
  
  open(vidObj);
  
  % Set up map data
  
  frame = zeros(resolution,2*resolution,3);
  
  RGCframe = []; %zeros(resolution,resolution,3);
  SCframe = zeros(resolution,resolution,3);
  
  % We want X to go N->T, and Y V->D
  
  [RGCnt,RGCdv] = meshgrid(linspace(1,0,resolution),linspace(1,0,resolution));
  [SCap,SCml] = meshgrid(linspace(0,1,resolution),linspace(0,1,resolution));
  
  fig = figure;
  for k = 1:numel(iterSequence)
    fprintf('Processing frame %d/%d (iter %d)\n', ...
            k, numel(iterSequence),iterSequence(k).curStep)

    [RGCcol,SCcol] = calculateSCcolour(iterSequence(k));

    frame = interpolateFrame(iterSequence(k),RGCcol,SCcol);
    imshow(frame);
    set(gca,'ydir','normal')
    
    % Add labels to retinal image
    textCol = [1 1 1]*0.95; %[209 27 141]/255;
    coordScale = resolution/500;
    
    text(resolution*0.5-15*coordScale,30*coordScale,'V', ...
         'fontsize',30*coordScale,'color',textCol);
    text(resolution*0.5-15*coordScale,resolution-25*coordScale,'D', ...
         'fontsize',30*coordScale,'color',textCol);
    text(15*coordScale,resolution*0.5,'T', ...
         'fontsize',30*coordScale,'color',textCol);
    text(resolution-40*coordScale,resolution*0.5,'N', ...
         'fontsize',30*coordScale,'color',textCol);

    text(resolution*1.5-15*coordScale,30*coordScale,'M', ...
         'fontsize',30*coordScale,'color',textCol);
    text(resolution*1.5-15*coordScale,resolution-25*coordScale,'L', ...
         'fontsize',30*coordScale,'color',textCol);
    text(resolution+40*coordScale,resolution*0.5,'A',...
         'fontsize',30*coordScale,'color',textCol);
    text(resolution*2-60*coordScale,resolution*0.5,'P',...
         'fontsize',30*coordScale,'color',textCol);
    
    text(resolution-70*coordScale,resolution-30*coordScale, ...
         sprintf('Iter = %d',iterSequence(k).curStep/iterSequence(k).nSC), ...
         'color', [1 1 1], 'fontsize',30*coordScale)
    
    % Mark a selected neuron
    switch(numel(SCidx))
      case 0
        % Do nothing
      case 1
        % Mark with white colour
        markSCneuron(k,SCidx,[1 1 1]);
      case 2
        markSCneuron(k,SCidx(1),[1 0 0]);
        markSCneuron(k,SCidx(2),[0 1 0]);
      otherwise
        for i = 1:numel(SCidx)
          markSCneuron(k,SCidx(i),[1 1 1]);
        end
    end

    % Mark a selected neuron
    switch(numel(RGCidx))
      case 0
        % Do nothing
      case 1
        % Mark with white colour
        markRGCneuron(k,RGCidx,[1 1 1]);
      case 2
        markRGCneuron(k,RGCidx(1),[1 0 0]);
        markRGCneuron(k,RGCidx(2),[0 1 0]);
      otherwise
        for i = 1:numel(RGCidx)
          markRGCneuron(k,RGCidx(i),[1 1 1]);
        end
    end
    
    drawnow
    % Print iteration
    finFrame = getframe(gcf);
    
    % Ugly hack, we need to find the interesting part of the data
    try
      img = finFrame.cdata;
      r = img(:,:,1);
      [y,x] = find(r == 0,1);
      finFrame.cdata = img(y:y+resolution-1,x:x+2*resolution-1,:);
      writeVideo(vidObj,finFrame);
    catch e
      getReport(e)
      keyboard
    end
  end

  fprintf('Writing movie to %s\n', movName)
  close(vidObj);
  
  toc
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function [RGCcol, SCcol] = calculateSCcolour(retObj)

    RGCb = (retObj.RGCnt - min(retObj.RGCnt)) ...
           / (max(retObj.RGCnt) - min(retObj.RGCnt));
    RGCg = (retObj.RGCdv - min(retObj.RGCdv)) ...
           / (max(retObj.RGCdv) - min(retObj.RGCdv));
  
    RGCcol = [zeros(size(RGCb)) RGCg RGCb];
    SCcol = zeros(retObj.nSC,3);
  
    for i = 1:retObj.nSC
      
      idx = retObj.presynapticConnections(1:retObj.numPresynapticConnections(i),i);
      w = retObj.presynapticWeight(1:retObj.numPresynapticConnections(i),i);
    
      for j = 1:numel(idx)
        SCcol(i,:) = SCcol(i,:) + RGCcol(idx(j),:)*w(j);
      end
    
      SCcol(i,:) = SCcol(i,:) / sum(w);
    
    end
    
    idx = find(isnan(SCcol(:,1)));
    SCcol(idx,:) = 0;
  
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function frame = interpolateFrame(retObj,RGCcol,SCcol)
   
    if(isempty(RGCframe))
      for i = 1:3
        RGCframe(:,:,i) = griddata(retObj.RGCnt,retObj.RGCdv,RGCcol(:,i), ...
                                   RGCnt,RGCdv);
      end
      
    end

    frame = zeros(resolution,2*resolution,3);

    for i = 1:3
      SCframe(:,:,i) = griddata(retObj.SCap,retObj.SCml,SCcol(:,i), ...
                                SCap,SCml);     
      
      frame(:,:,i) = [RGCframe(:,:,i),SCframe(:,:,i)];
    end
    
  end
  
  function markRGCneuron(iter, RGCidx, colour)
  
    iterSequence(iter).convertConnectionTables('pre2post');
    hold on
    
    [retX,retY] = ...
        coordToPixelPositionRetina(iterSequence(iter).RGCnt(RGCidx), ...
                                   iterSequence(iter).RGCdv(RGCidx));
    
    plot(retX,retY,'.','markersize',20,'color',colour);
    
    nCon = iterSequence(iter).numPostsynapticConnections(RGCidx);
    scIdx = iterSequence(iter).postsynapticConnections(1:nCon,RGCidx);
    SCW = iterSequence(iter).postsynapticWeight(1:nCon,RGCidx);
    
    [scX,scY] = coordToPixelPositionSC(iterSequence(iter).SCap(scIdx), ...
                                       iterSequence(iter) ...
                                       .SCml(scIdx));
    
    for i = 1:numel(scIdx)
      plot(scX(i),scY(i),'o','markersize',SCW(i)+3,'color',colour);
    end
    
    hold off
    
  end
    
  function markSCneuron(iter, SCidx, colour)
  
    hold on
    [scX,scY] = coordToPixelPositionSC(iterSequence(iter).SCap(SCidx), ...
                                       iterSequence(iter).SCml(SCidx));
    plot(scX,scY,'.','markersize',20,'color',colour);
    
    nCon = iterSequence(iter).numPresynapticConnections(SCidx);
    retIdx = iterSequence(iter).presynapticConnections(1:nCon,SCidx);
    retW = iterSequence(iter).presynapticWeight(1:nCon,SCidx);
    
    [retX,retY] = coordToPixelPositionRetina(iterSequence(iter).RGCnt(retIdx), ...
                                             iterSequence(iter).RGCdv(retIdx));
    
    for i = 1:numel(retIdx)
      plot(retX(i),retY(i),'o','markersize',retW(i)+3,'color',colour);
    end
    
    hold off
    
  end
    
  function [X,Y] = coordToPixelPositionRetina(nt,dv)
    
    % nt is Y (1st index)
    % dv is X (2nd index)
  
    % We got the X and Y axis flipped
    X = ceil((1-nt)*resolution);
    Y = ceil((1-dv)*resolution);
      
  end

  function [X,Y] = coordToPixelPositionSC(ap,ml)
    
    % ap is Y (1st index)
    % ml is X (2nd index)
    
    X = ceil(ap*resolution) + resolution;
    Y = ceil(ml*resolution);
    
  end
    
end