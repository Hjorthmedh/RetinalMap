% Plots one or two injections
% You need to specify SCidx and RGCidx, you can skip, or leave
% either of the two others empty. 

function [fig,sp] = plotVirtualInjection(obj,SCidx,RGCidx,SCidx2,RGCidx2)

    RGCmarker = '.';
    SCmarker = '.';
    RGC2marker = '.';
    SC2marker = '.';
    markerSize = 12;
    color1 = [1 0 0];
    color2 = [15 147 24]/255;
    otherColor = [1 1 1]*0.75; %0.8;
    
    if(~exist('SCidx2'))
      SCidx2 = [];
    end
    
    if(~exist('RGCidx2'))
      RGCidx2 = [];
    end

    fig = figure();
    
    sp(1) = subplot(1,2,1);
    idxFreeRGC = setdiff(1:obj.nRGC,union(RGCidx,RGCidx2));
    plot(obj.RGCnt(idxFreeRGC),obj.RGCdv(idxFreeRGC),...
         '.','color', otherColor);
    hold on
    plot(obj.RGCnt(RGCidx),obj.RGCdv(RGCidx),RGCmarker, ...
         'color',color1,'markersize',markerSize);
    plot(obj.RGCnt(RGCidx2),obj.RGCdv(RGCidx2),RGC2marker,...
         'color',color2,'markersize',markerSize);
    hold off
    xlabel('Temporal - Nasal')
    ylabel('Ventral - Dorsal')
    set(gca, 'xdir','reverse','ydir','reverse');
    axis equal
    axis([0 max(max(obj.RGCnt),1) 0 max(max(obj.RGCdv),1)])
    box off
    
    sp(2) = subplot(1,2,2);
    idxFreeSC = setdiff(1:obj.nSC,union(SCidx,SCidx2));
    plot(obj.SCap(idxFreeSC),obj.SCml(idxFreeSC), ...
         '.', 'color', otherColor);
    hold on
    plot(obj.SCap(SCidx),obj.SCml(SCidx),SCmarker, ...
         'color',color1,'markersize',markerSize);
    plot(obj.SCap(SCidx2),obj.SCml(SCidx2),SC2marker,...
         'color',color2,'markersize',markerSize);
    hold off
    xlabel('Anterior - Posterior')
    ylabel('Medial - Lateral')
    axis equal
    axis([0 max(max(obj.SCap),1) 0 max(max(obj.SCml),1)])
    box off
    

end