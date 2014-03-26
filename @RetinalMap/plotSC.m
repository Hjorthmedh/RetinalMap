function plotSC(obj,noAxisFlag)

  assert(strcmp(obj.eyeType,'disk'))

  if(~exist('noAxisFlag'))
    noAxisFlag = 0;
  end
  
  % Calculate the SC colours
  RGCb = (obj.RGCnt - min(obj.RGCnt)) / (max(obj.RGCnt) - min(obj.RGCnt));
  RGCg = (obj.RGCdv - min(obj.RGCdv)) / (max(obj.RGCdv) - min(obj.RGCdv));
  
  RGCcol = [zeros(size(RGCb)) RGCg RGCb];
  SCcol = zeros(obj.nSC,3);
  
  for i = 1:obj.nSC
    idx = obj.presynapticConnections(1:obj.numPresynapticConnections(i),i);
    w = obj.presynapticWeight(1:obj.numPresynapticConnections(i),i);
    
    for j = 1:numel(idx)
      SCcol(i,:) = SCcol(i,:) + RGCcol(idx(j),:)*w(j);
    end
    
    SCcol(i,:) = SCcol(i,:) / sum(w);
    
  end
  
  % Plot the SC neurons
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
  

    if(~noAxisFlag)
      xlabel('Anterior - Posterior','fontsize',16)
      ylabel('Medial - Lateral','fontsize',16)
      set(gca,'fontsize',16)
     
      axis equal
      axis([0 1 0 1])
    else
      axis equal
      axis([0 1 0 1])
      axis off
    end
    
    box off

end
