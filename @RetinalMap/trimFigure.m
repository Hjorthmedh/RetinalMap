function trimFigure(obj,fig)

  % This function trims off unecessary whitespace around a plot

  ch = get(fig,'children');
  
  if(isempty(ch) | nnz(ismember(get(ch,'type'),{'axes'})) == 0)
    % No figure to trim!
    return;
  end
  
  [minx,maxx,miny,maxy,...
   minxMargin,maxxMargin,minyMargin,maxyMargin] = getMinMax(fig,'normalized');
  
  % Update all positions accordingly

  for i = 1:numel(ch)
    switch(get(ch(i),'type'))
      case 'axes'
        pos = get(ch(i),'position');

        % We want x1 y1 x2 y2, not x1 y1 w h
        pos([3 4]) = pos([1 2]) + pos([3 4]);
        
        x = [1 3]; y = [2 4];
        pos(x) = ((pos(x) - minx) / (maxx-minx)) ...
                 * (1 - minxMargin - maxxMargin) + minxMargin;
        pos(y) = ((pos(y) - miny) / (maxy-miny)) ...
                 * (1 - minyMargin - maxyMargin) + minyMargin;

        pos([3 4]) = pos([3 4]) - pos([1 2]);

        set(ch(i),'position',pos);
      case {'uimenu','uitoolbar','uicontextmenu'}
        % Ignore these
      otherwise
        fprintf('Unhandled type: %s\n', get(ch(i),'type'))
    end
  end

  set(fig,'units','centimeters');

  [minx,maxx,miny,maxy,...
   minxMargin,maxxMargin,minyMargin,maxyMargin] = getMinMax(fig,'centimeters');
  
  wh = [maxx-minx+minxMargin+maxxMargin, ...
        maxy-miny+minyMargin+maxyMargin];

  set(fig,'paperunits','centimeters');
  set(fig,'papersize', wh);
  set(fig,'paperpositionmode','manual');
  set(fig,'paperposition',[0, 0, wh]);
  
  set(fig,'units','normalized');
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function [minx,maxx,miny,maxy,...
            minxMargin,maxxMargin,minyMargin,maxyMargin] = getMinMax(fig,units)

    minxMargin = 0;
    minyMargin = 0;
    maxxMargin = 0;
    maxyMargin = 0;
    
    minx = inf;
    miny = inf;
    maxx = -inf;
    maxy = -inf;
    
    % Find min and max extent
    for i = 1:numel(ch)
      switch(get(ch(i),'type'))
        case 'axes'
    
          %oldUnits = get(ch(i),'units');
          
          set(ch(i),'units',units);
          pos = get(ch(i),'position');
          ti = get(ch(i),'tightinset');
          
          %set(ch(i),'units',oldUnits);

          
          % We want x1 x2 y1 y2, not x1 y1 w h
          pos([3 4]) = pos([1 2]) + pos([3 4]);
        
          if(pos(1) < minx)
            minx = pos(1);
            minxMargin = ti(1);
          elseif(pos(1) == minx)
            minxMargin = max(minxMargin,ti(1));
          end
          
          if(pos(2) < miny)
            miny = pos(2);
            minyMargin = ti(2);
          elseif(pos(2) == miny)
            minyMargin = max(minyMargin,ti(2));
          end
          
          if(pos(3) > maxx)
            maxx = pos(3);
            maxxMargin = ti(3);
          elseif(pos(3) == maxx)
            maxxMargin = max(maxxMargin,ti(3));
          end
          
          if(pos(4) > maxy)
            maxy = pos(4);
            maxyMargin = ti(4);
          elseif(pos(4) == maxy)
            maxyMargin = max(maxyMargin,ti(4));
          end
          
        otherwise
          % Do nothing
      end
    end
  end
  
end