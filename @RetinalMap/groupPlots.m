% This function takes existing figures and places them in subplots
% to better use the space

function fig = groupPlots(obj,subplotY,subplotX,figHandleList,singleFlag,...
                          fontsize,markersize,markerSizeDot,linewidth, ...
                          subplotPositionIdx)

  if(~exist('singleFlag'))
    singleFlag = ones(size(figHandleList));
  end

  if(numel(singleFlag) == 1)
    singleFlag = singleFlag*ones(size(figHandleList));
  end
  
  if(obj.plotFigures == 0)
    disp('groupPlots: plotFigures = 0, hiding figures!')
    visFlag = 'off';
  else
    visFlag = 'on';
  end

  fig = figure('visible',visFlag);
  set(fig,'papertype','A4');

  set(fig,'paperunits','centimeters');
  set(fig,'paperposition',[0 0 21 29.7])
  %set(fig,'paperposition',[0.6345 0.6345 20.3046 28.431])
  set(fig,'units','centimeters');
  set(fig,'position',[0 0 21 29.7]);  
  set(fig,'units','normalized');

  
  
  % Trim the figures
  for i = 1:numel(figHandleList)
    obj.trimFigure(figHandleList(i));
  end
  
  % First make a list of all axes that should be in each subplot
  axesList = {};
  expandAxis = [];
  
  for i = 1:numel(figHandleList)
    ch = allchild(figHandleList(i));
    idx = find(ismember(get(ch,'type'),{'axes'}));

    if(singleFlag(i))
      axesList{end+1} = ch(idx);
      expandAxis(end+1) = false;
    else
      for j = numel(idx):-1:1
        axesList{end+1} = ch(idx(j));
        
        if(numel(idx) > 1)
          % Only need to resize it if it is not single plot
          expandAxis(end+1) = true;
        else
          expandAxis(end+1) = false;
        end
      end
    end
  end
  
  if(numel(axesList) > subplotX*subplotY)
    fprintf('%d x %d subplot can not fit %d plots\n', ...
            subplotY, subplotX, numel(axesList))
    subplotY = ceil(numel(axesList)/subplotX);
  end

  for i = 1:numel(axesList)
    if(exist('subplotPositionIdx'))
      sp = subplot(subplotY,subplotX,subplotPositionIdx(i));
    else
      % Default behaviour
      sp = subplot(subplotY,subplotX,i);
    end

    spPos = get(sp,'position');    
    delete(sp);
      
    for j = numel(axesList{i}):-1:1
      set(axesList{i}(j),'units','normalized')
      oldPos = get(axesList{i}(j),'position');
      newAxes = copyobj(axesList{i}(j),fig);
      
      try
        if(expandAxis(i))
          newPos = spPos;
        else
          newPos = [spPos(1) + spPos(3)*oldPos(1); 
                    spPos(2) + spPos(4)*oldPos(2);
                    spPos(3)*oldPos(3);
                    spPos(4)*oldPos(4)];
        end

        if(exist('fontsize'))
          set(newAxes,'fontsize',fontsize)
          ch = allchild(newAxes);
          idx = find(ismember(get(ch,'type'),{'text'}));
          set(ch(idx),'fontsize',fontsize)
        end      
      
        set(newAxes,'position',newPos);      
      
      catch e
        getReport(e)
        keyboard
      end

      set(newAxes,'Xdir', get(axesList{i}(j),'Xdir'));
      set(newAxes,'Ydir', get(axesList{i}(j),'Ydir'));
      set(newAxes,'visible', get(axesList{i}(j),'visible'));
      
      % colormap(get(get(axesList{i}(j),'parent'),'colormap'))
      
      
      if(exist('markersize'))
        idx = find(ismember(get(ch,'type'),{'line'}));
        set(ch(idx),'markersize',markersize);
        
        try
          idx2 = idx(find(ismember(get(ch(idx),'marker'),'.')));
          set(ch(idx2),'markersize',markerSizeDot);
        catch e
          getReport(e)
          keyboard
        end
      end

      if(exist('linewidth'))
        idx = find(ismember(get(ch,'type'),{'line'}));
        set(ch(idx),'linewidth',linewidth)
      end
      
      % !!! Get colourmap from old figure
      % colormap('hot')
      
    end
            
  end
  
  % obj.trimFigure(fig);
  
  % keyboard
  
end
