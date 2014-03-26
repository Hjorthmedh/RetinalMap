% Code by David Sterratt
%
% Normally SCcol and RGCcol are not specified by the user

function fig = plotGradients2D(obj, invertSCephrinA, forwardOnly, markIsl2)

  silentFlag = 1;

  if (exist('forwardOnly') && ~forwardOnly) 
    forwardOnly = false;
  else 
    forwardOnly = true;
  end
  if (exist('markIsl2') && markIsl2)
    markIsl2 = 1;
  else 
    markIsl2 = 0;
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % This allows the user to invert the SC ephrinA gradient

  SCephrinA = obj.SCephrinA;
  if (exist('invertSCephrinA') && invertSCephrinA)
    SCephrinA = 1./obj.SCephrinA;
  end
  
  minA = min([obj.RGCEphA ; SCephrinA]);
  maxA = max([obj.RGCEphA ; SCephrinA]);
  minB = min([obj.RGCEphB ; obj.SCephrinB]);
  maxB = max([obj.RGCEphB ; obj.SCephrinB]);
  
  % Calculate colour of each SC neuron
  SCA = (SCephrinA - minA)/(maxA - minA);
  SCB = (obj.SCephrinB - minB)/(maxB - minB); 

  %
  RGCA = (obj.RGCEphA - minA)/(maxA - minA);
  RGCB = (obj.RGCEphB - minB)/(maxB - minB);

  
  if(obj.plotFigures == 0)
    disp('plotMapReverse: plotFigures = 0, hiding figures!')
    visFlag = 'off';
  else
    visFlag = 'on';
  end
  
  fig = figure('visible',visFlag);
  
  if (forwardOnly) 
    subplot(2,2,1); 
  else 
    subplot(2,4,1);
  end 
  axesRetinaA = gca();
  plotRetinalMap(getColour(RGCA));
  title('Retinal EphA')
  
  if (forwardOnly) 
    subplot(2,2,2); 
  else 
    subplot(2,4,2);
  end 
  axesSCA = gca();
  plotSCmap(getColour(SCA));
  title('SC ephrin-A')

  if (forwardOnly) 
    subplot(2,2,3);
  else 
    subplot(2,4,5);
  end 
  axesRetinaB = gca();
  plotRetinalMap(getColour(RGCB))
  title('Retinal EphB')
  
  if (forwardOnly) 
    subplot(2,2,4); 
  else 
    subplot(2,4,6);
  end 
  axesSCB = gca();
  plotSCmap(getColour(SCB));
  title('SC ephrin-B')

  if (~forwardOnly) 
    subplot(2,4,3);
    axesRetinaA = gca();
    plotRetinalMap(getColour(obj.RGCephrinA));
    title('Retinal ephrin-A')
    
    subplot(2,4,4);
    axesSCA = gca();
    plotSCmap(getColour(obj.SCEphA));
    title('SC EphA')
    
    subplot(2,4,7);
    axesRetinaB = gca();
    plotRetinalMap(getColour(obj.RGCephrinB))
    title('Retinal ephrin-B')
    
    subplot(2,4,8);
    axesSCB = gca();
    plotSCmap(getColour(obj.SCEphB));
    title('SC EphB')
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function plotRetinalMap(RGCcol)
    
    switch(obj.eyeType)
      case 'sphere'
        hold off
        polar(linspace(0,2*pi,100),obj.maxTheta*ones(1,100),'r-')
        hold on
  
        for i = 1:obj.nRGC
          r = obj.RGCtheta(i);
          v = obj.RGCphi(i);
        
          if(any(isnan(RGCcol(i,:))))
            plot(r*cos(v),r*sin(v), 'ko');
          else
            plot(r*cos(v),r*sin(v), '.', 'color', RGCcol(i,:), ...
                 'markersize',12);
          end
        end

        title('Remember to add axes labels')
      
        axis equal         
      
      case 'disk'
        hold off
        if(markIsl2 & ~isempty(obj.Isl2PositiveRGC))
          plot(obj.RGCnt(obj.Isl2PositiveRGC), ...
               obj.RGCdv(obj.Isl2PositiveRGC), ...
               'o','color',[1 0 0]*0.8, ...
               'markersize', 6);
          hold on
        end
        for i = 1:obj.nRGC
          if(any(isnan(RGCcol(i,:))))
            plot(obj.RGCnt(i), obj.RGCdv(i), 'ko');
          else
            plot(obj.RGCnt(i), obj.RGCdv(i),...
                 '.', 'color', RGCcol(i,:), ...
                 'markersize', 12);
          end
          hold on        
        end

        % We want T-N on X axis, and V-D on Y-axis, need to reverse them
        set(gca, 'xdir','reverse','ydir','reverse');        
        
        ylabel('Ventral - Dorsal','fontsize',16)
        xlabel('Temporal - Nasal','fontsize',16)
        set(gca,'fontsize',16)
    
        axis equal    
        axis([0 max(max(obj.RGCnt),1) 0 max(max(obj.RGCdv),1)])
        box off
            
      otherwise
        fprintf('plotMapReverse: Unknown eye type %s\n', obj.eyeType)
        keyboard  
    end
    
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function plotSCmap(SCcol)

    hold off
    for i = 1:numel(obj.SCml)
      if(obj.numPresynapticConnections(i) > 0)
        plot(obj.SCap(i),obj.SCml(i),'.', ...
             'color',SCcol(i,:), ...
             'markersize',12);
      else
        plot(obj.SCap(i),obj.SCml(i),'ko');        
      end
      hold on
    end
  

    xlabel('Anterior - Posterior','fontsize',16)
    ylabel('Medial - Lateral','fontsize',16)
    set(gca,'fontsize',16)
  
    if(0)
      % Mark the axons
      
      for i = 1:numel(obj.RGCml)
        plot([0 1],obj.RGCml(i)*[1 1],'k-')
      end
    end
    
    axis equal
    axis([0 max(max(obj.SCap),1) 0 max(max(obj.SCml),1)])
    box off
  
  end
  
end

function RGB = getColour(X)
  X = (X - min(X))/(max(X) - min(X));
  RGB = [0.1 + 0.9*X, 0.5 + 0.5*X, ... 
         repmat(0.7, size(X, 1), 1)];
end
