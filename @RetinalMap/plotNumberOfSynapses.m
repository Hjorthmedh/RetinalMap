function fig = plotNumberOfSynapses(obj,plotType,oldFig)

  if(obj.plotFigures == 0)
    disp('plotNumberOfSynapses: plotFigures = 0, hiding figures!')
    visFlag = 'off';
  else
    visFlag = 'on';
  end
    
  if(exist('oldFig') & ~isempty(oldFig))
    % This is currently only supported for the cumdist plots
    assert(strcmp(plotType,'cumdist'));
    reuseFigFlag = true;
  else
    reuseFigFlag = false;
  end
  
  switch(plotType)
  
    case 'map'
    
      maxSyn = max(max(obj.totalWeightRGC),max(obj.totalWeightSC));

      maxSynR = max(obj.totalWeightRGC);
      maxSynSC = max(obj.totalWeightSC);
  
      c = colormap('gray');
    
      fig = figure('visible',visFlag);
  
      subplot(1,2,1)
      hold on
      for i = 1:obj.nRGC
        cR = interp1(linspace(0,maxSynR,size(c,1)),c(:,1), ...
                     obj.totalWeightRGC(i));
        cG = interp1(linspace(0,maxSynR,size(c,1)),c(:,2), ...
                     obj.totalWeightRGC(i));
        cB = interp1(linspace(0,maxSynR,size(c,1)),c(:,3), ...
                     obj.totalWeightRGC(i));
        plot(obj.RGCdv(i),obj.RGCnt(i), ...
             '.','color', [cR cG cB], 'markersize', 12)
      end
      axis equal
      axis([0 max(max(obj.RGCnt),1) 0 max(max(obj.RGCdv),1)])
      box off
      
      colormap('gray');
      cb = colorbar;
      yVal = 0:10:maxSynR;
      ytick = interp1(linspace(0,maxSynR,size(c,1)),1:size(c,1),yVal);
      yTickLabel = {};
    
      for i = 1:numel(yVal)
        yTickLabel{i} = num2str(yVal(i));
      end
  
      set(cb,'ytick',ytick,'yticklabel',yTickLabel) 
      set(gca,'xdir','reverse','ydir','reverse')
      xlabel('Temporal-Nasal')
      ylabel('Ventral-Dorsal')
      box off
      
      subplot(1,2,2)
      hold on
      for i = 1:obj.nSC
        try
          cR = interp1(linspace(0,maxSynSC,size(c,1)),c(:,1), ...
                       obj.totalWeightSC(i));
          cG = interp1(linspace(0,maxSynSC,size(c,1)),c(:,2), ...
                       obj.totalWeightSC(i));
          cB = interp1(linspace(0,maxSynSC,size(c,1)),c(:,3), ...
                       obj.totalWeightSC(i));
          plot(obj.SCap(i),obj.SCml(i), ...
               '.','color', [cR cG cB], 'markersize', 12)
        catch e
          getReport(e)
          keyboard
        end
      end
      
      cb = colorbar;
      yVal = 0:2:maxSynSC;
      ytick = interp1(linspace(0,maxSynSC,size(c,1)),1:size(c,1),yVal);
      yTickLabel = {};
      for i = 1:numel(yVal)
        yTickLabel{i} = num2str(yVal(i));
      end
  
      set(cb,'ytick',ytick,'yticklabel',yTickLabel) 
  
      axis equal
      axis([0 max(max(obj.SCap),1) 0 max(max(obj.SCml),1)])
      xlabel('Anterior-Posterior')
      ylabel('Medial-Lateral')
      
      box off
  
      fName = sprintf('%s/%s-numberOfSynapses-map.pdf', ...
                      obj.figurePath,obj.simName);
  
      switch(obj.plotFigures)
        case {0,1}
          saveas(fig,fName,'pdf');
        case 2
          obj.addParametersToFigure(strrep(fName,'.pdf','.eps'));
      end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    case 'hist'
  
      obj.convertConnectionTables('pre2post');
  
      fig = figure('visible',visFlag);
      subplot(2,2,1)
  
      nTotPreSynCon = sum(obj.presynapticWeight);
      nTotPostSynCon = sum(obj.postsynapticWeight);
  
      plot(obj.RGCnt,nTotPostSynCon,'o','color',[1 1 1]*0.7)
      medianPlot(obj.RGCnt,nTotPostSynCon);
      xlabel('Nasal-Temporal')
      ylabel('Number of synapses')
      title('Retina')
      box off
  
      subplot(2,2,3)
      plot(obj.RGCdv,nTotPostSynCon,'o','color',[1 1 1]*0.7)
      medianPlot(obj.RGCdv,nTotPostSynCon);
      xlabel('Dorsal-Ventral')
      ylabel('Number of synapses')
      box off

      subplot(2,2,2)
      plot(obj.SCap,nTotPreSynCon,'o','color',[1 1 1]*0.7)
      medianPlot(obj.SCap,nTotPreSynCon);
      xlabel('Anterior-Posterior')
      ylabel('Number of synapses')  
      title('Superior colliculus')
      box off
  
      subplot(2,2,4)
      plot(obj.SCml,nTotPreSynCon,'o','color',[1 1 1]*0.7)
      medianPlot(obj.SCml,nTotPreSynCon);
      xlabel('Medial-Lateral')
      ylabel('Number of synapses')  
      box off
      
      fName = sprintf('%s/%s-numberOfSynapses-plot.pdf', ...
                      obj.figurePath,obj.simName);
      
      switch(obj.plotFigures)
        case {0,1}
          saveas(fig,fName,'pdf');
        case 2
          obj.addParametersToFigure(strrep(fName,'.pdf','.eps'));
      end
      

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      
    case 'cumdist'

      % Data extracted from Tripplet et al 2011, figure 3D
      
      WTlowerAP = load('math5/Triplett2011-figs-3D-lower-WT.csv');      
      WTupperAP = load('math5/Triplett2011-figs-3D-upper-WT.csv');
      Math5lowerAP = load('math5/Triplett2011-figs-3D-lower-Math5.csv');
      Math5upperAP = load('math5/Triplett2011-figs-3D-upper-Math5.csv');

      WTlowerML = load('math5/Triplett2011-figs-3E-lower-WT.csv');      
      WTupperML = load('math5/Triplett2011-figs-3E-upper-WT.csv');
      Math5lowerML = load('math5/Triplett2011-figs-3E-lower-Math5.csv');
      Math5upperML = load('math5/Triplett2011-figs-3E-upper-Math5.csv');
      
      if(reuseFigFlag)
        fig = oldFig;
        figure(fig);
        set(fig,'visible',visFlag);
        sp = get(fig,'children');
      else
        fig = figure('visible',visFlag);
      end
      
      [~,idxAP] = sort(obj.SCap);
      [~,idxML] = sort(obj.SCml);
  
      try
        cumWeightSCap = cumsum(obj.totalWeightSC(idxAP));
        cumWeightSCml = cumsum(obj.totalWeightSC(idxML));
      catch e
        getReport(e)
        keyboard
      end
      
      if(reuseFigFlag)
        set(fig,'currentaxes',sp(2));
      else
        subplot(2,1,1)

        if(strcmpi(obj.phenotype,'WT'))
          p = patch([WTupperAP(:,1);WTlowerAP(end:-1:1,1)], ...
                    [WTupperAP(:,2);WTlowerAP(end:-1:1,2)], [1 0.4 0.4]*0.7);
          set(p,'edgecolor',get(p,'facecolor'));
        elseif(strcmpi(obj.phenotype,'Math5'))
          p = patch([Math5upperAP(:,1);Math5lowerAP(end:-1:1,1)], ...
                    [Math5upperAP(:,2);Math5lowerAP(end:-1:1,2)], [1 0.4 0.4]*0.7);
          set(p,'edgecolor',get(p,'facecolor'));
        end
      end
      
      hold on
      
      s = stairs(obj.SCap(idxAP),cumWeightSCap/max(cumWeightSCap));
      set(s,'color',[0 0 0]);
      hold on
      %plot([min(obj.SCap) max(obj.SCap)],[0 1],'k-')
      % xlabel('Anterior-Posterior','fontsize',24)
      set(gca,'xtick',[0.1 0.9],'xticklabel',{'A','P'});
      set(gca,'ticklength', [0 0]);

      ylabel('Cumulative #synapses (normalised)','fontsize',24)
      set(gca,'fontsize',30)
      box off
      
      title('Includes all SC neurons')

      if(reuseFigFlag)
        set(fig,'currentaxes',sp(1));
      else
        subplot(2,1,2)

        if(strcmpi(obj.phenotype,'WT'))
          p = patch([WTupperML(:,1);WTlowerML(end:-1:1,1)], ...
                    [WTupperML(:,2);WTlowerML(end:-1:1,2)], [1 0.4 0.4]*0.7);
          set(p,'edgecolor',get(p,'facecolor'));
        elseif(strcmpi(obj.phenotype,'Math5'))
          p = patch([Math5upperML(:,1);Math5lowerML(end:-1:1,1)], ...
                    [Math5upperML(:,2);Math5lowerML(end:-1:1,2)], [1 0.4 0.4]*0.7);
          set(p,'edgecolor',get(p,'facecolor'));        
        end
      end
      
      hold on
      
      s = stairs(obj.SCml(idxML),cumWeightSCml/max(cumWeightSCml));
      set(s,'color',[0 0 0]);  
      hold on
      %plot([min(obj.SCml) max(obj.SCml)],[0 1],'k-')
      % xlabel('Medial-Lateral','fontsize',24)
      set(gca,'xtick',[0.1 0.9],'xticklabel',{'M','L'});
      set(gca,'ticklength', [0 0]);
      
      ylabel('Cumulative #synapses (normalised)','fontsize',24)
      set(gca,'fontsize',30)
      
      box off
      
      fName = sprintf('%s/%s-cumulative-normalisedNumberOfSynapses-SC.pdf', ...
                      obj.figurePath,obj.simName);
      
      switch(obj.plotFigures)
        case {0,1}
          saveas(fig,fName,'pdf');
        case 2
          obj.addParametersToFigure(strrep(fName,'.pdf','.eps'));
      end
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    case 'time'
  
      % Plot the number of synapses per RGC and SC neuron nas a
      % function of time
      
      fig = figure('visible',visFlag);
      
      skipTables = 1;
      s = obj.loadIterSequence(skipTables);
      
      iter = cat(1,s.curStep);
      medianNumSynapsesRGC = zeros(numel(s),1);
      p5NumSynapsesRGC = zeros(numel(s),1);
      p95NumSynapsesRGC = zeros(numel(s),1);
      
      for i = 1:numel(s)
        medianNumSynapsesRGC(i) = median(s(i).totalWeightRGC);
        p5NumSynapsesRGC(i) = prctile(s(i).totalWeightRGC,5);
        p95NumSynapsesRGC(i) = prctile(s(i).totalWeightRGC,95);        
      end

      a = area(repmat(iter,1,2)/obj.nSC, ...
           [p5NumSynapsesRGC, ...
            p95NumSynapsesRGC-p5NumSynapsesRGC], ...
          'facecolor', 0.6*[1 1 1], ...
          'edgecolor', 0.6*[1 1 1]);
  
      delete(a(1)); 
      hold on      
      
      plot(iter/obj.nSC,medianNumSynapsesRGC,'k-','linewidth',2)
      xlabel('Iteration','fontsize',20)
      ylabel('Number of synapses per RGC','fontsize',20)
      set(gca,'fontsize',20)
      
      fName = sprintf('%s/%s-numberOfSynapses-per-iteration.pdf', ...
                      obj.figurePath,obj.simName);
       
      switch(obj.plotFigures)
        case {0,1}
          saveas(fig,fName,'pdf');
        case 2
          obj.addParametersToFigure(strrep(fName,'.pdf','.eps'));
      end
      
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  function [xBin, yMed, y5, y95] = medianPlot(x,y)
    
    binSize = 0.1;
    xBin = 0:binSize:1;
    
    yMed = zeros(size(xBin));
    y5 = zeros(size(xBin));
    y95 = zeros(size(xBin));   
    
    for i = 1:numel(xBin)
      idx = find(xBin(i)-binSize/2 <= x ...
                 & x < xBin(i) + binSize/2);
      
      yMed(i) = median(y(idx));
      y5(i) = prctile(y(idx),5); 
      y95(i) = prctile(y(idx),95); 
      
    end
    
    hold on
    plot(xBin,yMed,'r-','linewidth',3);
    plot(xBin,y5,'r--','linewidth',1);
    plot(xBin,y95,'r--','linewidth',1);    
    
    hold off
  end
  
end