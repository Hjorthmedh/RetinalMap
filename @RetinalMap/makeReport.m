function report = makeReport(obj, plotFigures, DWgrid)

  if(~exist('plotFigures'))
    plotFigures = 1;
  else
    r.plotFigures = plotFigures;
  end

  if(~exist('DWgrid'))
    disp('Not using DW grid')
    DWgrid = false;
  end

  fontSize = 4;  
  
  % Return some info to the caller
  report.simName = obj.simName;     
 
  allFigures = [];
  singlePlotFlag = [];
 
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Running DWs code, it does not play nicely with existing
  %% figures, will always overwrite figure 2,3 and 6.
      
  if(DWgrid)
      
    fakeID = ceil(rand*1000000);
    fprintf('Using fakeID = %d\n', fakeID)
    [DWresults,DWfigHandle] = willshawGrid(obj,fakeID);
    
    if(isempty(DWresults))
      disp('DW grid analysis failed, trying to recover.')
      % Something went wrong, try to recover, exclude DW grid from analysis
      DWgrid = false;  
    end
    % Add report figures further down
  else
    DWresults = [];
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  summaryInit();
  
  figPar = obj.plotParameters(fontSize);
  addToSummaryFig(figPar,1);

  fig = obj.plotGradients1D();
  addToSummaryFig(fig(1),1);
  addToSummaryFig(fig(2),1);      
      
  fig = obj.plotMapForward(1);
  addToSummaryFig(fig,1);
      
  fig = obj.plotMapImage();
  addToSummaryFig(fig,1);

  % Find largest uncrossed submap
  [report.grid,fig] = obj.makeProjectionGrid();
  addToSummaryFig(fig,1);  
  
  fig = newFig();
  obj.plotAxisMap('AP',[],[],1,plotFigures)

  APcollapsePoint = obj.locateCollapsePoint();
  if(0 < APcollapsePoint & APcollapsePoint < 1)
    hold on        
    a = axis();
    plot(APcollapsePoint*[1 1]*(a(2)-a(1))+a(1),a(3:4),'w--')
    hold off
    title(sprintf('Collapse point %.2f', APcollapsePoint))
  end
  addToSummaryFig(fig,1);
  
  report.APcollapsePoint = APcollapsePoint;
  
  fig = newFig();
  obj.plotAxisMap('ML',[],[],1,plotFigures)
  addToSummaryFig(fig,1);
  
  fig = newFig();
  obj.plotAxisMap('EphA',[],[],1,plotFigures)      
  addToSummaryFig(fig,1);      
  
  
  fig = obj.virtualInjectionMappingExperiment('NT')
  addToSummaryFig(fig,1);

  if(0)
    disp('Plotting map fuzzyness - calculating dFuzz can take some time')
    [dFuzzNT,kFuzzAP,cFuzz95AP,fig] = obj.plotMapFuzzyness('NT');
    addToSummaryFig(fig,1);
  
  
    [dFuzzDV,kFuzzML,cFuzz95ML,fig] = obj.plotMapFuzzyness('DV');
    addToSummaryFig(fig,1);
  
    report.dFuzzNT = dFuzzNT;
    report.kFuzzAP = kFuzzAP;
    report.cFuzz95AP = cFuzz95AP;
  
    report.dFuzzDV = dFuzzDV;
    report.kFuzzML = kFuzzML;
    report.cFuzz95ML = cFuzz95ML;
  end
    
  % Termination zone size
  fig = obj.plotTerminationZoneSize();
  if(~isempty(fig))
    addToSummaryFig(fig(1),1);
    addToSummaryFig(fig(2),1);  
  end
  
  % Virtual retrograde injections
  [report.kSegregationAP, fig] = ...
      obj.plotVirtualExperimentSegregation('AP');
  addToSummaryFig(fig,1);
  
  [report.kSegregationML,fig] = ...
      obj.plotVirtualExperimentSegregation('ML');
  addToSummaryFig(fig,1);
  
  [report.kSegregationSeqAP,report.iterSegregationAP,fig] = ...
      obj.plotVirtualExperimentSegregationSequence('AP');
  addToSummaryFig(fig,1);
  
  [report.kSegregationSeqML,report.iterSegregationML,fig] = ...
      obj.plotVirtualExperimentSegregationSequence('ML');
  addToSummaryFig(fig,1);  
  
  fig = obj.plotVirtualExperimentVector('AP');
  addToSummaryFig(fig,1);
  
  SCidx = obj.findSCneuron(0.5,0.5,0.028/2);
  [RGCidx,weight] = obj.findPresynapticRGC(SCidx);

  plotContours = true;
  
  [areaContour, ~, fig] = obj.contourAnalysisKDE(RGCidx,[0.05 0.25 0.5 0.75 0.95],plotContours);  
  addToSummaryFig(fig,1);


  report.areaContour = areaContour;
  
  % Find superposed projection in SC
  [report.superposedProj, ...
   report.superposedProjMajorAng, ...
   report.superposedProjMajorMinorRatio, ...
   fig] = obj.superposedProjection();
  addToSummaryFig(fig,1);
  
  %%% Additional check for Isl2
  
  if(1)
    if(~isempty(obj.Isl2PositiveRGC))
      % Make two additional plots for Isl2 positive and negative
      % RGC
      
      [report.superposedProjIsl2Pos, ...
       report.superposedProjMajorAngIsl2Pos, ...
       report.superposedProjMajorMinorRatioIsl2Pos, ...
       fig] = obj.superposedProjection(obj.Isl2PositiveRGC);
      addToSummaryFig(fig,1);
      
      [report.superposedProjIsl2Neg, ...
       report.superposedProjMajorAngIsl2Neg, ...
       report.superposedProjMajorMinorRatioIsl2Neg, ...
       fig] = obj.superposedProjection(setdiff(1:obj.nRGC,obj.Isl2PositiveRGC));
      addToSummaryFig(fig,1);
      
    else
      % No Isl2-EphA3, lets leave these empty
      
      report.superposedProjIsl2Pos = [];
      report.superposedProjMajorAngIsl2Pos = [];
      report.superposedProjMajorMinorRatioIsl2Pos = [];
      
      report.superposedProjIsl2Neg = [];
      report.superposedProjMajorAngIsl2Neg = [];
      report.superposedProjMajorMinorRatioIsl2Neg = [];
      
    end
  end
  
  %%% End Isl2
  
  % !!! Do postage stamp plots, to make sure that there are not
  % two populations of projections, for say Isl2.

  % Check if nasal or temporal injections give ectopics
  try
    [nasalEctopic,fig] = obj.virtualInjectionEctopicExperiment(0.2,0.5);
    addToSummaryFig(fig,1);
    [temporalEctopic,fig] = obj.virtualInjectionEctopicExperiment(0.8,0.5);
    addToSummaryFig(fig,1);      
  catch e
    getReport(e)
    keyboard
  end
  
  fprintf('Nasal ectopic: %d, Temporal ectopic: %d\n', ...
          nasalEctopic, temporalEctopic)
  
  
  report.nasalEctopic = nasalEctopic;
  report.temporalEctopic = temporalEctopic;

  %% Number of synapses
  
  fig = obj.plotNumberOfSynapses('cumdist'); %'map','hist','cumdist'
  addToSummaryFig(fig,1);
  
  %%%%%%%
  
  report.SC95coverage = obj.calculateSCcoverage(0.95);
  
  
  % Store the entire simulation in the report structure... for now
  obj.plotFigures = 1;
  report.RetinalMap = obj;
  
  if(DWgrid)
    try
      report.DWgrid = DWresults;
      if(0)
        % Include DW figures from main summary page?
        for iFig = 1:numel(DWfigHandle)
          addToSummaryFig(DWfigHandle(iFig),1);
        end
      end
    catch e
      getReport(e)
      keyboard
    end
  else
    report.DWgrid = [];
  end
  
  fig = obj.plotQuantificationPanel(report,5);
  addToSummaryFig(fig,1);
  
  % Add a report sheet also with quantified data??
  
  % Make a summary plot
  fig = obj.groupPlots(3,3,allFigures,singlePlotFlag,5,1,3,1);
  
  fName = sprintf('%s/%s-one-page-summary.eps',obj.figurePath,obj.simName);
  fprintf('Printing to file %s\n', fName)
  print(fig,'-dpsc2',fName)
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Make a second result figure, with bigger figures

  if(DWgrid)

    summaryInit();

    for iFig = 1:numel(DWfigHandle)
      addToSummaryFig(DWfigHandle(iFig),0);
    end
    
    % Add parameters
    addToSummaryFig(figPar,1);    
    
    fig2 = obj.groupPlots(6,3,allFigures,singlePlotFlag,...
                          6,5,1,1,[1:6 7:8 10:11 13:14 16:17 18]);

    fprintf('Appending to file %s\n', fName)
    obj.saveMultiPagePS(fName,fig2,1,'eps');
  
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function fig = newFig()
    if(plotFigures)
      fig = figure();
    else
      fig = figure('visible','off');
    end
  end
    
  function addToSummaryFig(fig,singlePlot)
    
    if(isempty(fig))
      disp('No figure to add, creating a blank figure')
      fig = newFig();
    end
    
    try
      for i = 1:numel(fig)
        allFigures(end+1) = fig(i);
        
        if(numel(singlePlot) > 1)
          singlePlotFlag(end+1) = singlePlot(i);          
        else
          singlePlotFlag(end+1) = singlePlot;
        end
      end
    catch e
      getReport(e)
      keyboard
    end
  end
  
  function summaryInit()
    allFigures = [];
    singlePlotFlag = [];
  end
  
end
