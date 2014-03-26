% The function plotMapForward colours the SC neurons depending on
% the centroid of the presynaptic neurons, and plotMapReverse
% colours the retinal neurons depending on the centroid of the post
% synaptic projection.
%
% This function instead plots all the termination centroids for the
% RGC, and colours them with the colour of to the RGC neuron.
%
% This function is useful to illustrate the Math5 mapping, since
% when using plotMapForward a single synapse on a SC neuron makes
% it coloured, even though the centroid for the incoming RGC might
% be far away.

function plotAxonTerminationCentroids(obj,onlyPlotSCflag)

  if(~exist('onlyPlotSCflag'))
    onlyPlotSCflag = 0;
  end
  
  assert(strcmp(obj.eyeType,'disk'))

  % Colour the RGC
  RGCb = (obj.RGCnt - min(obj.RGCnt)) / (max(obj.RGCnt) - min(obj.RGCnt));
  RGCg = (obj.RGCdv - min(obj.RGCdv)) / (max(obj.RGCdv) - min(obj.RGCdv));

  RGCcol = [zeros(size(RGCb)) RGCg RGCb];

  nCon = zeros(obj.nRGC,1);
  RGCcentAP = zeros(obj.nRGC,1);
  RGCcentML = zeros(obj.nRGC,1);
  
  % Find the SC projection centroids of all the RGC
  for i = 1:obj.nSC
    for j = 1:obj.numPresynapticConnections(i)
      idx = obj.presynapticConnections(j,i);
      w = obj.presynapticWeight(j,i);
      RGCcentAP(idx) = RGCcentAP(idx) + w*obj.SCap(i);
      RGCcentML(idx) = RGCcentML(idx) + w*obj.SCml(i);
      nCon(idx) = nCon(idx) + w;
    end
  end
  
  % Internal check, these should match
  assert(nnz(nCon-obj.totalWeightRGC) == 0);
  
  RGCcentAP = RGCcentAP ./ nCon;
  RGCcentML = RGCcentML ./ nCon;
  
  if(~onlyPlotSCflag)
    fig = figure;
    
    subplot(1,2,1);
    retinaAxes = gca();
    % Plot the RGC
    for i = 1:obj.nRGC
      plot(obj.RGCnt(i),obj.RGCdv(i),'.', ...
           'markersize',10,'color',RGCcol(i,:));
      hold on
    end
    hold off
  
    set(gca,'xdir','reverse','ydir','reverse');
    xlabel('Temporal-Nasal','fontsize',18)
    ylabel('Ventral-Dorsal','fontsize',18)
    set(gca,'fontsize',16)
    axis equal
    axis([0 1 0 1]);
    box off


    subplot(1,2,2)
    scAxes = gca();
  
  end

  % Plot the outline of SC
  idx = convhull(obj.SCap,obj.SCml);
  plot(obj.SCap(idx),obj.SCml(idx),'k--')
  hold on
  
  % Plot the centroids in the RGC colour

  for i = 1:obj.nRGC
    plot(RGCcentAP(i),RGCcentML(i),'.', ...
         'color', RGCcol(i,:), 'markersize',10);
  end
  hold off
  xlabel('Anterior-Posterior','fontsize',18)
  ylabel('Medial-Lateral','fontsize',18)
  set(gca,'fontsize',16)
  axis equal
  axis([0 1 0 1]);
  box off

  if(~onlyPlotSCflag)
    fName = sprintf('%s/%s-termination-centroids-iter-%d.pdf', ...
                    obj.figurePath,obj.simName,obj.curStep/obj.nSC);
  
    if(~exist(obj.figurePath))
      mkdir(obj.figurePath)
    end
  
    saveas(gcf,fName,'pdf');

    set(fig,'currentaxes',retinaAxes)
    set(gca,'color',0.2*[1 1 1]);  
    set(fig,'currentaxes',scAxes)
  end
  
  set(gca,'color',0.2*[1 1 1]);  
  
end