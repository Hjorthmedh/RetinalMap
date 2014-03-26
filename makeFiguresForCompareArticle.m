% This code will get integrated into ComparePhenotype code, but
% right now I am running simulations, so want to make sure that
% code does not get an error introduced, and stop running half way through.
%
% makeFiguresForCompareArticle('WT',1001:1010,'koulakov')
% makeFiguresForCompareArticle('Isl2homozygous',1001:1010,'koulakov')
% makeFiguresForCompareArticle('Isl2heterozygous',1001:1010,'koulakov')
% makeFiguresForCompareArticle('TKO',1001:1010,'koulakov')
% makeFiguresForCompareArticle('ephrinA2mm',1001:1010,'koulakov')
% makeFiguresForCompareArticle('ephrinA5mm',1001:1010,'koulakov')
% makeFiguresForCompareArticle('Math5',1001:1010,'koulakov')
%
% makeFiguresForCompareArticle('WT',2001:2010,'koulakov')
% makeFiguresForCompareArticle('Isl2homozygous',2001:2010,'koulakov')
% makeFiguresForCompareArticle('Isl2heterozygous',2001:2010,'koulakov')
% makeFiguresForCompareArticle('TKO',2001:2010,'koulakov')
% makeFiguresForCompareArticle('ephrinA2mm',2001:2010,'koulakov')
% makeFiguresForCompareArticle('ephrinA5mm',2001:2010,'koulakov')
% makeFiguresForCompareArticle('Math5',2001:2010,'koulakov')

% makeFiguresForCompareArticle('WT',[4,6:12,100,101],'Markerinduction')
% makeFiguresForCompareArticle('Isl2homozygous',[4,6:12,100,101],'Markerinduction')
% makeFiguresForCompareArticle('Isl2heterozygous',[4,6:12,100,101],'Markerinduction')
% makeFiguresForCompareArticle('TKO',[4,6:12,100,101],'Markerinduction')
% makeFiguresForCompareArticle('Math5',[4,6:12,100,101],'Markerinduction')



function report = makeFiguresForCompareArticle(phenotype,expRange,model,plotFigures)
  
  if(~exist('model') | isempty(model))
    model = 'koulakov';
  end
  
  fontSize = 4;  
  markerSize = 1;
  markerSizeDot = 3;
  lineWidth = 1;
  
  dataDir = 'SAVE/ComparePhenotype';
  
  for i = 1:numel(expRange)
    simName{i} = sprintf('%s/ComparePhenotype-%s-%s-rep-%d.mat', ...
                         dataDir, model, phenotype, expRange(i));
  end
  
  disp('Loading data')
  
  %% Load the data
  for i = 1:numel(simName)
    memoryFlag = 1;
    r(i) = RetinalMap();
    r(i).loadState(simName{i},memoryFlag);
    
    % Debugging purpose
    if(exist('plotFigures') & ~isempty(plotFigures))
      r(i).plotFigures = plotFigures;
    end
  end
  
  if(~exist(r(1).figurePath))
    mkdir(r(1).figurePath);
  end
  
  expStr = sprintf('-%d',expRange);
  summaryFigName = sprintf('FIGS/ComparePhenotype/ComparePhenotype-%s-%s-summary%s.eps',...
                           model, phenotype, expStr);

  reportName = sprintf('FIGS/ComparePhenotype/ComparePhenotype-%s-%s-report%s.mat',...
                           model, phenotype, expStr);
  
  
  switch(lower(phenotype))
    case 'wt'

      %plotList = {'params','gradients1D','map','picmap','gridJH',...
      %            'NT','segregation','contours','ectopic','synapses','summary'};
      
      
      plotList = {'gradients1D','map','gridJH', ...
                  'NThist','DVhist','NTinj', ...
                  'segregation','contours','synapses','ectopic',...
                  'DWgrid','params','summary'};
  
    case {'isl2homozygous','isl2heterozygous'}

      plotList = {'gradients1D','map','gridJH',...
                  'NThist','DVhist','NTinj','segregation','DWgrid','params','summary'};

      % plotList = {'NThist', 'DWgrid'};
      
    case 'math5'

      plotList = {'gradients1D','map','gridJH',...
                  'NThist','DVhist','segregation','synapses', 'wholeeye', ...
                  'DWgrid','params','summary'};
      
      
    case 'tko'

      plotList = {'gradients1D','map','gridJH',...
                  'NThist','DVhist','ectopic','DWgrid','params','summary'};
      
      
      % plotList = {'NThist', 'DVhist', 'DWgrid'};
      
    case {'ephrina2mm','ephrina5mm'}
      
      plotList = {'gradients1D','map','gridJH',...
                  'NThist','DVhist','ectopic','DWgrid','params','summary'};
      
      
    case {'isl2hom', 'isl2het', ...
                'isl2hetepha4mm', 'isl2homepha4mm', ...
                'isl2hetepha4pm', 'isl2homepha4pm', ...
                'isl2hetepha5mm', 'isl2homepha5mm', ...
                'isl2hetepha5pm', 'isl2homepha5pm'}
      
      plotList = {'NThist'};
      
    otherwise
      
      fprintf('Unknown phenotype: %s\n', phenotype)
      keyboard
      
  end    
    
  % In case we just want to do one analysis for all phenotypes...
  %plotList = {'DWgrid'}
  % plotList = {'segregation'}
  % plotList = {'NThist'}
  % plotList = {'DVhist'}  
  % plotList = {'synapses'}
  % plotList = {'wholeeye'}
  
  summaryInit();
  
  if(r(1).plotFigures)
    visFlag = 'on';
  else
    visFlag = 'off';
  end
    
  report = [];
  
  report.phenotype = phenotype;
  report.model = model;
  
  switch(model)
    case 'WhiteCow'
      weightScale = 0.075;
      alpha = 0.05;
    case 'Markerinduction'
      weightScale = 10;
      alpha = 0.02;
    otherwise
      weightScale = 1;
      alpha = 0.05;
  end
  
  for plotIdx = 1:numel(plotList)
    
    switch(plotList{plotIdx})
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      case 'params'

        % Just displaying the first experiments parameters
        fig = r(1).plotParameters(fontSize);
        addToSummaryFig(fig,1);

        % Check that all parameters are equal, all that are
        % different will be displayed in the command window.
        for i = 2:numel(r)
          r(1).paramDiff(r(i));
        end

        
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
      case 'gradients1D'
        
        % Just showing one example of gradients

        fig = r(1).plotGradients1D();
        addToSummaryFig(fig(1),1);
        addToSummaryFig(fig(2),1);      
          
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      case 'map'
        
        for i = 1:numel(r)
          %fig = r(i).plotMapForward(1);
          switch(lower(phenotype))
            case 'tko'
              % We want to more easilly see the "islands" of order
              fig = r(i).plotMapImage('images/grid-red-green-blue-yellow.png');              
            otherwise
              fig = r(i).plotMapImage('images/grid-red-green-blue-yellow-black.png');
          end
          
          %addToSummaryFig(fig,0);
          addToSummaryFig(fig,1);
        end
        
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      case 'picmap'
        
        for i = 1:numel(r)
          fig = r(i).plotMapImage();
          addToSummaryFig(fig,1);
        end
        
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      case 'gridJH'
        
        for i = 1:numel(r)
          [report.grid(i),fig] = r(i).makeProjectionGrid();
          % addToSummaryFig(fig,0);  
          addToSummaryFig(fig,1);            
        end
          
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      case 'NThist'

        plotType = 'hist'; % 'alpha';
                           % weightScale and alpha ignored if using
                           % option 'hist'
        
        showAll = false; % Only show middle third
        reuseFigFlag = false;
        useAlpha = true;
        saveFigFlag = true;
        
        switch(lower(phenotype))
          case {'isl2heterozygous', ...
                'isl2hom', 'isl2het', ...
                'isl2hetepha4mm', 'isl2homepha4mm', ...
                'isl2hetepha4pm', 'isl2homepha4pm', ...
                'isl2hetepha5mm', 'isl2homepha5mm', ...
                'isl2hetepha5pm', 'isl2homepha5pm'}
        
            APcollapsePoint = zeros(numel(r),1);
        
            for i = 1:numel(r)
              fig = r(i).plotAxisProjection('NT',showAll,reuseFigFlag, ...
                                              plotType,saveFigFlag, ...
                                              weightScale,alpha);
              
              APcollapsePoint(i) = r(i).locateCollapsePoint();
              if(0 < APcollapsePoint(i) & APcollapsePoint(i) < 1)
                hold on        
                a = axis();
                
                plot(APcollapsePoint(i)*[1 1]*(a(2)-a(1))+a(1),a(3:4),'k--')                  
                hold off
                title(sprintf('Collapse point %.0f %%', 100*APcollapsePoint(i)))
              
                % Overwrite the previous saved
                fName = sprintf('%s/%s-NT-projection.pdf', ...
                                r(i).figurePath, r(i).simName);
                
                print(gcf,'-dpdf',fName,'-painters','-r1200');
              
              end

              
              addToSummaryFig(fig,1);
            end
        
            
            report.APcollapsePoint = APcollapsePoint;
            
          otherwise

            for i = 1:numel(r)
              fig = r(i).plotAxisProjection('NT',showAll,reuseFigFlag, ...
                                            plotType,saveFigFlag, ...
                                            weightScale,alpha);
              
              addToSummaryFig(fig,1);
            end
            
            report.APcollapsePoint = [];
        end
          
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
      case 'DVhist'

        plotType = 'hist';

        showAll = false; % Only show middle third
        reuseFigFlag = false;
        useAlpha = true;
        saveFigFlag = true;
        
        for i=1:numel(r)
          fig = r(i).plotAxisProjection('DV',showAll,reuseFigFlag,...
                                        plotType,saveFigFlag,weightScale,alpha);
          addToSummaryFig(fig,1);
        end
        
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      case 'NTinj'

        for i = 1:numel(r)
          fig = r(i).virtualInjectionMappingExperiment('NT');
          addToSummaryFig(fig,1);
        end
        
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      case 'segregation'

        % Displays the segregation curves, AP direction
        % Dashed line, WT P60 comparison
        % Solid lines, segregation fits to virtual injections
        
        % Get experimental fits
        scaleDist = true;        
        visSegFlag = false;
        
        rHelper = RetinalMap();
        [fitAge, kWT, kB2KO, kWTJK, kB2KOJK] = rHelper.fitSegregationData(scaleDist,visSegFlag);

        nInj = 200;
        injectionType = 'MLfix';

        xRange = linspace(0,1,100);
        logistic = inline('1 - 0.5./(1+(x/k(1)).^k(2))','k','x');
        
        % Get segregation fits
        for i = 1:numel(r)
          [~,~,report.kFitSegregation(:,i)] = ...
              r(i).virtualInjectionSegregationExperiment(injectionType,nInj);
        end

        % Get Jack-knife error bars
        yAll60 = zeros(numel(xRange),size(kWTJK{end},2));
        for i = 1:size(kWTJK{end},2)
          yAll60(:,i) = logistic(kWTJK{end}(:,i),xRange);
        end
        yMax60 = max(yAll60,[],2);
        yMin60 = min(yAll60,[],2);

        yAll22 = zeros(numel(xRange),size(kWTJK{end-1},2));
        for i = 1:size(kWTJK{end-1},2)
          yAll22(:,i) = logistic(kWTJK{end-1}(:,i),xRange);
        end
        yMax22 = max(yAll22,[],2);
        yMin22 = min(yAll22,[],2);

        
        
        % Add check for visible flag in r(1)
        fig = figure('visible',visFlag);
        hold on
        a = area(transpose(xRange), [yMin60, yMax60-yMin60], ...
                 'facecolor', 0.6*[1 1 1] + 0.4*[1 0 0], ...
                 'edgecolor', 0.6*[1 1 1] + 0.4*[1 0 0]);
        delete(a(1));
         
        a = area(transpose(xRange), [yMin22, yMax22-yMin22], ...
                 'facecolor', 0.6*[1 1 1] + 0.4*[1 0 0], ...
                 'edgecolor', 0.6*[1 1 1] + 0.4*[1 0 0]);
        delete(a(1));
       
       
        switch(lower(phenotype))
          case 'wt'
            plot(xRange,logistic(kWT(end,:),xRange),'r-','linewidth',2);
            plot(xRange,logistic(kWT(end-1,:),xRange),'r--','linewidth',2);          
          otherwise
            plot(xRange,logistic(kWT(end,:),xRange),'b-','linewidth',2);            
            plot(xRange,logistic(kWT(end-1,:),xRange),'b--','linewidth',2);                    
        end
            
        yAll = zeros(numel(xRange),numel(r));
        
        for i = 1:numel(r)
          yAll(:,i) = logistic(report.kFitSegregation(:,i),xRange);
        end
        
        plot(xRange,yAll, '-','color',0.6*[1 1 1]);
        plot(xRange,median(yAll,2),'k-')
        
        xlabel('Normalised SC distance','fontsize',24)    
        ylabel('Segregation in retina','fontsize',24)
        set(gca,'fontsize',24)
        set(gca,'ytick',0.5:0.1:1)
        box off
        
        addToSummaryFig(fig,1);

        % Save a plot also
        expStr = sprintf('-%d',expRange);
        fName = sprintf('%s/ComparePhenotype-%s-%s-segregation-summary%s.pdf', ...
                        r(1).figurePath,model,phenotype,expStr);
        saveas(fig, fName, 'pdf');
        
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
      case 'contours'
        
        % Load experimental data from Dan
        refWT = importdata('contours/WT-P60.csv');
        
        meanContour = nanmean(refWT.data(:,3:7));
        stdContour = nanstd(refWT.data(:,3:7));
        
        percentile = [0.95 0.75 0.5 0.25 0.05];
        nRep = 20;
        contourPlotFlag = false;
        
        allContours = zeros(nRep,numel(percentile),numel(r));
        meanContours = zeros(numel(percentile),numel(r));
        
        for i = 1:numel(r)
          allContours(:,:,i) = ...
                      r(i).calculateTerminationZoneSizeInjection(nRep,percentile,contourPlotFlag);
        
          meanContours(:,i) = mean(allContours(:,:,i));
        end
        
        fig = figure('visible',visFlag);
        errorbar(100*percentile,meanContour,stdContour,'r--')
        hold on
        plot(100*repmat(transpose(percentile),1,numel(r)),meanContours,'ko')        

        xlabel('Percentage of retinal labeling','fontsize',24)
        ylabel('Fraction of retina labeled','fontsize',24)
        set(gca,'fontsize',24)
        box off  
        
        addToSummaryFig(fig,1);
        
        % Save a plot also
        expStr = sprintf('-%d',expRange);
        fName = sprintf('%s/ComparePhenotype-%s-%s-contours-summary%s.pdf', ...
                        r(1).figurePath,model,phenotype,expStr);
        saveas(fig, fName, 'pdf');

        
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
      case 'wholeeye'
        
        % Do whole eye labeling, how much is covered by 95% of label?
        
        for i = 1:numel(r)
          [areaCov(i)] = r(i).calculateSCcoverage(0.99);
        end
        
        report.areaCoverage = areaCov;
        

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      case 'ectopic'
      
        for i = 1:numel(r)
          [nasalEctopic(i),fig] = r(i).virtualInjectionEctopicExperiment(0.1,0.5);
          % addToSummaryFig(fig,0);
          addToSummaryFig(fig,1);    
        end

        for i = 1:numel(r)        
          [temporalEctopic(i),fig] = r(i).virtualInjectionEctopicExperiment(0.9,0.5);
          % addToSummaryFig(fig,0);      
          addToSummaryFig(fig,1);                
        end
        
        report.nasalEctopic = nasalEctopic;
        report.temporalEctopic = temporalEctopic;
        
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

      case 'synapses'
        
        % Math5 cumulative synapse plot
        
        fig = [];
        for i = 1:numel(r)
          fig = r(i).plotNumberOfSynapses('cumdist',fig); %'map','hist','cumdist'
        end
        
        addToSummaryFig(fig,1);
        
        % Save a plot also
        expStr = sprintf('-%d',expRange);
        fName = sprintf('%s/ComparePhenotype-%s-%s-synapse-cumulative-distribution-summary%s.pdf', ...
                        r(1).figurePath,model,phenotype,expStr);
        saveas(fig, fName, 'pdf');
         
         
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

      case 'summary'
      
        disp('summary: Not yet implemented')
        
        
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      case 'DWgrid'

        savedFig = zeros(1,numel(r));
        
        submapSizeNodes = zeros(numel(r),2)*NaN;
        submapSizeEdges = zeros(numel(r),2)*NaN;        
        
        for i = 1:numel(r)
          
          fakeID = ceil(rand*1000000);
          fprintf('Using fakeID = %d\n', fakeID)
          [DWresults,DWfigHandle,P] = willshawGrid(r(i),fakeID);
    
          if(isempty(DWresults))
            disp('DW grid analysis failed, trying to recover.')
          end
          
          % Each time we run the lattice analysis, we overwrite the
          % figures. Not good, so we got to copy them to a new
          % figure to preserve them. 
          %
          % Update: Not needed if we use DS plot 78.
          
          for j = 1:numel(DWfigHandle)
            savedFig(j,i) = figure('visible',visFlag);
            copyobj(get(DWfigHandle(j),'children'),savedFig(j,i));
          end

          submapSizeNodes(i,1) = 100*(1 - ...
                                        P(1).stats.FTOC.num_nodes_crossing ...
                                        / P(1).FTOC.numpoints);
          
          submapSizeEdges(i,1) = P(1).FTOC.percent_edges_in_subgraph;          

          if(numel(P) == 2)
            submapSizeNodes(i,2) = 100*(1 - ...
                                        P(2).stats.FTOC.num_nodes_crossing ...
                                        / P(2).FTOC.numpoints);
            
            submapSizeEdges(i,2) = P(2).FTOC.percent_edges_in_subgraph;          
          end
          
        end

        % Sort the figures and output them one type at a time
        for iFig = 1:size(savedFig,1)
          for i = 1:numel(r)
            addToSummaryFig(savedFig(iFig,i),1);
          end
        end
          
        report.submapSizeNodes = submapSizeNodes;
        report.submapSizeEdges = submapSizeEdges;
        
        
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        

      otherwise
        
        fprintf('Unknown plot requested: %s\n', plotList{plotIdx})
        keyboard
        
    end
    
  end

  % Save report
  save(reportName,'report');
  
  % Make a summary plot
  appendFlag = false;  
  nRows = 7; %8;
  nCols = 3; %4;
  nPlots = nRows*nCols; % Desired #figs per sheet

  figStart = 1;
  figEnd = 0;
  
  while(figEnd < numel(allFigures))
    nCtr = 0; % Figures on current sheet, as of now
    figStart = figEnd+1;

    if(singlePlotFlag(figEnd+1))
      n = 1;
    else
      n = numel(get(allFigures(figEnd+1),'children'));
    end
    
    while(nCtr+n <= nPlots & figEnd < numel(allFigures))
      nCtr = nCtr + n;
      figEnd = figEnd + 1;
      
      if(figEnd >= numel(allFigures))
        % We are done here
        continue
      end
      
      % We got to count the axes that this fix will add to the figure
      if(singlePlotFlag(figEnd+1))
        n = 1;
      else
        n = numel(get(allFigures(figEnd+1),'children'));
      end
    end

    figIdx = figStart:figEnd;
    fig = r(1).groupPlots(nRows,nCols,allFigures(figIdx),singlePlotFlag(figIdx),...
                          fontSize,markerSize,markerSizeDot,lineWidth);
    r(1).saveMultiPagePS(summaryFigName,fig,appendFlag,'eps');
    appendFlag = true;

    figEnd = figEnd + 1;
    
    % keyboard
    
  end
      
  if(~r(1).plotFigures)
    disp('Invisible figures, closing all.')
    close all
  end
    
  % Convert eps to pdf
  outFile = strrep(summaryFigName,'.eps','.pdf');
  sysCmd = sprintf('ps2pdf %s %s', summaryFigName,outFile);
  system(sysCmd);
  
  % SJE: Now crop the pdf to a tighter bounding box.
  sysCmd = sprintf('pdfcrop %s', outFile);
  system(sysCmd);  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function summaryInit()
    allFigures = [];
    singlePlotFlag = [];
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function fig = newFig()
    fig = figure('visible',visFlag);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
end