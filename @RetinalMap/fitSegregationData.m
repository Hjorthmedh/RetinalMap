function [age, kWT, kB2KO, kWTJK, kB2KOJK] = fitSegregationData(obj,scaleDist,visibleFlag)

  if(~exist('scaleDist') | isempty(scaleDist))
    scaleDist = false;
  end
  
  if(~exist('visibleFlag'))
    visibleFlag = 'on';
  else
    if(visibleFlag)
      visibleFlag = 'on';
    else
      visibleFlag = 'off';
    end
  end
  
  useLoess = false;
  
  nRep = 100;
  nGrid = 100;  
  
  distance = [];
  segregation = [];
  
  SCsize = [];
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % We use Dans parameters as starting points for nlinfit
  [WTfit,B2KOfit,SCsize] = obj.getDanSegregationFits();
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % These files are for the AP axis
  % Note that distScale uses AP axis size also (see function below)
  % WTfile = 'segregation/segregation-WT-C57BL6J.csv';
  % B2KOfile = 'segregation/segregation-B2KO.csv';
  
  % Latest AP files from Dan Lyngholm, 2012-12-12
  WTfile = 'segregation/NN-C57BL6J-AP.csv';
  B2KOfile = 'segregation/NN-B2KO-AP.csv';

  
  WTdata = importdata(WTfile);
  B2KOdata = importdata(B2KOfile);
  
  WT.age = WTdata.data(:,1);
  WT.dist = WTdata.data(:,2);
  WT.distScaled = scaleSCdist(WT.dist,WT.age);
  WT.NNred = WTdata.data(:,3)/100;
  WT.NNgreen = WTdata.data(:,4)/100;
  WT.NNmean = WTdata.data(:,5)/100;

  B2KO.age = B2KOdata.data(:,1);
  B2KO.dist = B2KOdata.data(:,2);
  B2KO.distScaled = scaleSCdist(B2KO.dist,B2KO.age);
  B2KO.NNred = B2KOdata.data(:,3)/100;
  B2KO.NNgreen = B2KOdata.data(:,4)/100;
  B2KO.NNmean = B2KOdata.data(:,5)/100;
  
  age = union(WT.age,B2KO.age);
  col = winter(numel(age));
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if(scaleDist)
    xMax = max(max(WT.distScaled),max(B2KO.distScaled));
  else
    xMax = max(max(WT.dist),max(B2KO.dist));
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  nPar = 2;
  
  kWT = NaN*zeros(numel(age),nPar);
  kB2KO = NaN*zeros(numel(age),nPar);

  nWT = zeros(numel(age),1);
  nB2KO = zeros(numel(age),1);
  
  kWTJK = {};
  kB2KOJK = {};
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % This does jack-knifing
  
  if(1)
   for i = 1:numel(age)

      figure('visible',visibleFlag)
      idxWT = find(WT.age == age(i));
      idxB2KO = find(B2KO.age == age(i));
      
      idxFitWT = find(WTfit.age == age(i));
      idxFitB2KO = find(B2KOfit.age == age(i));
      
      if(scaleDist)
        distanceWT = WT.distScaled(idxWT);      
        distanceB2KO = B2KO.distScaled(idxB2KO);   
        k0WT = WTfit.kAPscaled(idxFitWT,:);
        k0B2KO = B2KOfit.kAPscaled(idxFitB2KO,:);
      else
        distanceWT = WT.dist(idxWT);
        distanceB2KO = B2KO.dist(idxB2KO);
        k0WT = WTfit.kAP(idxFitWT,:);
        k0B2KO = B2KOfit.kAP(idxFitB2KO,:);
      end  
      
      segregationWT = WT.NNmean(idxWT);
      segregationB2KO = B2KO.NNmean(idxB2KO);
      
      if(numel(idxWT) > 3)
        
        if(useLoess)
          [xRangeWT,yMedianWT,yAllWT] = obj.loessJackKnife(distanceWT,segregationWT);
        else
          [kWT(i,:),kWTJK{i},xRangeWT,yMedianWT,yAllWT] = ...
              obj.fitJackKnife(distanceWT,segregationWT,@logistic,k0WT);
        end
      else
        xRangeWT = [];
        yMedianWT = [];
        yAllWT = [];
      end    
      
      if(numel(idxB2KO) > 3)
        if(useLoess)
          [xRangeB2KO,yMedianB2KO,yAllB2KO] = obj.loessJackKnife(distanceB2KO,segregationB2KO);          
        else
          [kB2KO(i,:),kB2KOJK{i},xRangeB2KO,yMedianB2KO,yAllB2KO] = ...
              obj.fitJackKnife(distanceB2KO,segregationB2KO,@ ...
                               logistic,k0B2KO);
        end
      else
        xRangeB2KO = [];
        yMedianB2KO = [];
        yAllB2KO = [];      
      end    
          
      hold on
      
      if(~isempty(xRangeWT))
        a = area(transpose(xRangeWT), ...
                 [min(yAllWT,[],2), ...
                  max(yAllWT,[],2)-min(yAllWT,[],2)], ...
                 'facecolor', 0.6*[1 1 1], ...
                 'edgecolor', 0.6*[1 1 1]);
        delete(a(1));
        plot(xRangeWT,yAllWT,'-','color',[1 1 1]*0.6);
        p(1) = plot(xRangeWT,yMedianWT,'k-','linewidth',2);
      end
        
      if(~isempty(xRangeB2KO))
        a = area(transpose(xRangeB2KO), ...
                 [min(yAllB2KO,[],2), ...
                  max(yAllB2KO,[],2)-min(yAllB2KO,[],2)], ...
                 'facecolor', 0.4*[1 0 0] + 0.6*[1 1 1], ...
                 'edgecolor', 0.4*[1 0 0] + 0.6*[1 1 1]);
        delete(a(1));        
        plot(xRangeB2KO,yAllB2KO, ...
             '-','color',[1 1 1]*0.6 + 0.4*[1 0 0]);
        p(2) = plot(xRangeB2KO,yMedianB2KO,'r-','linewidth',2);
      end
      
      plot(distanceWT,segregationWT,'k.', ...
           distanceB2KO,segregationB2KO,'r.', ...
           'markersize',20);
      
      if(obj.curStep > 0)
        % Lets add the simulation data also
        [distanceModel,segregationModel,kModel] = ...
            obj.virtualInjectionSegregationExperiment('MLfix',100);
        
        x = linspace(0,1);
        p(3) = plot(x,logistic(kModel,x),'-b');
        plot(distanceModel,segregationModel,'b.');
        legend(p,'WT','B2KO','model','location','southeast');

      else
        legend(p,'WT','B2KO','location','southeast');
      end
      

   
      xlabel('SC distance','fontsize',20)
      ylabel('Segregation in retina','fontsize',20)
      set(gca,'fontsize',20)
      
      if(~scaleDist)
        if(age(i) <= 6)
          axMax = 2000;
        else
          axMax = 1000;
        end
        axis([0 axMax 0.4 1])
      end
      
      % axis([0 xMax 0 1])
      
      if(obj.curStep > 0)
        title(sprintf('Age P%d (iter %d)', age(i),obj.curStep/obj.nSC))

        fName = sprintf('%s/%s-Dan-segregation-data-P%d-jackknife.pdf', ...
                  obj.figurePath, obj.simName, age(i));

      else
        title(sprintf('Age P%d', age(i)))
        fName = sprintf('FIGS/Dan-segregation-data-P%d-jackknife.pdf',age(i));
      end
      saveas(gcf,fName,'pdf')
   
   end

  end
  
  % Cluster and plot parameter space for different ages
  plotParamSpace();
  
  plotParamSpaceALT();  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % This does bootstrapping
  if(0)
    for i = 1:numel(age)

      figure('visible',visibleFlag)
      idxWT = find(WT.age == age(i));
      idxB2KO = find(B2KO.age == age(i));
      
      if(scaleDist)
        distanceWT = WT.distScaled(idxWT);      
        distanceB2KO = B2KO.distScaled(idxB2KO);   
        k0 = [rand(1),2*rand(1)];
      else
        distanceWT = WT.dist(idxWT);
        distanceB2KO = B2KO.dist(idxB2KO);
        k0 = [1600*rand(1),2*rand(1)];    
      end  
      
      segregationWT = WT.NNmean(idxWT);
      segregationB2KO = B2KO.NNmean(idxB2KO);
      
      if(numel(idxWT) > 3)
        [xRangeWT,yMedianWT,yLowWT,yHighWT,yAllWT] = ...
            obj.fitErrorEstimateN(distanceWT,segregationWT,@logistic,k0);
      else
        xRangeWT = [];
        yMedianWT = [];
        yLowWT = [];
        yHighWT = [];
        yAllWT = [];
      end    
      
      if(numel(idxB2KO) > 3)
        [xRangeB2KO,yMedianB2KO,yLowB2KO,yHighB2KO,yAllB2KO] = ...
            obj.fitErrorEstimateN(distanceB2KO,segregationB2KO,@logistic,k0);
      else
        xRangeB2KO = [];
        yMedianB2KO = [];
        yLowB2KO = [];
        yHighB2KO = [];
        yAllB2KO = [];      
      end    
      
      
      hold on
      if(~isempty(xRangeWT))
        a = area(repmat(xRangeWT,1,2), ...
                 [yLowWT,yHighWT-yLowWT], ...
                 'facecolor', 0.6*[1 1 1], ...
                 'edgecolor', 0.6*[1 1 1]);
        
        pWT = plot(xRangeWT,yMedianWT,'k-','linewidth',2);
        
        delete(a(1)); 
        alpha(get(a(2),'children'),0.5);
      end
      
      if(~isempty(xRangeB2KO))
        a = area(repmat(xRangeB2KO,1,2), ...
                 [yLowB2KO,yHighB2KO-yLowB2KO], ...
                 'facecolor', 0.4*[1 0 0] + 0.6*[1 1 1], ...
                 'edgecolor', 0.4*[1 0 0] + 0.6*[1 1 1]);
        
        pB2KO = plot(xRangeB2KO,yMedianB2KO,'r-','linewidth',2);
        
        delete(a(1)); 
        alpha(get(a(2),'children'),0.5);    
        
      end
      
      plot(distanceWT,segregationWT,'k.', ...
           distanceB2KO,segregationB2KO,'r.', ...
           'markersize',15);
      
      legend([pWT,pB2KO],'WT','B2KO');
      xlabel('SC distance','fontsize',20)
      ylabel('Segregation in retina','fontsize',20)
      set(gca,'fontsize',15)
      
      title(sprintf('Age P%d', age(i)))
      % axis([0 xMax 0 1])
      
      fName = sprintf('FIGS/Dan-segregation-data-P%d.pdf',age(i));
      saveas(gcf,fName,'pdf')
      
    end
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function y = logistic(k,x)
    % y = -0.5./(1+(x/inflection).^slope)+1;
    y = 1 - 0.5./(1+(x/k(1)).^k(2));
    
    % If k(1) is negative, we can get imaginary solutions
    % If k(2) is negative, x=0 --> y = 1.
    if(any(k < 0))
      y = ones(size(x))*1e5;
    end
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function x = inverseLogistic(k,y)
    if(numel(y) == 1)
      x = k(:,1).*((0.5-y)./(y-1)).^(1./k(:,2));    
    else
      x = k(1)*((0.5-y)./(y-1)).^(1/k(2));
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function dScaled = scaleSCdist(dRaw,age)
    % We want the largest axis to be between 0 and 1, so scale by AP size 
  
    SCage = SCsize(:,1);
    SCAP = SCsize(:,2);
    SCML = SCsize(:,3);
    
    scaleFactor = interp1(SCage,SCAP,age);
    
    dScaled = dRaw ./ scaleFactor;
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function plotParamSpace()
    
    col = [228, 26, 28;
           55, 126, 184;
           77, 175, 74;
           152, 78, 163;
           255, 127, 0;
           255, 255, 51;
           166, 86, 40;
           247, 129, 191;
           153, 153, 153]/255;
           
    col = [255, 255, 217;
           237, 248, 217; 
           199, 233, 180; 
           127, 205, 187; 
           65, 182, 196; 
           29, 145, 192; 
           34, 94, 168; 
           37, 52, 148; 
           8, 29, 88]/255;
    
    
    nWT = zeros(numel(age),1);
    kWTall = [];
    kWTage = [];
    nB2KO = zeros(numel(age),1);
    kB2KOall = [];
    
    figure('visible',visibleFlag)
    subplot(2,1,1)
    for i = 1:numel(age)
      p(i) = loglog(kWTJK{i}(1,:),kWTJK{i}(2,:),'.','color',col(i,:));
      hold on
      % plot(kB2KOJK{i}(1,:),kB2KOJK{i}(2,:),'.','color',col(i,:));      
      pLeg{i} = sprintf('P%d',age(i));
      
      nWT(i) = size(kWTJK{i},2);
      kWTall = [kWTall,kWTJK{i}];
      kWTage = [kWTage,age(i)*ones(1,nWT(i))];
    end

    xlabel('Inflection point')
    ylabel('Slope')
    legend(p,pLeg,'location','west');
    title('Wildtype')
    clear p pLeg
    
    subplot(2,1,2)
    nClustWT = 3;
    
    [idxWT, centroidsWT] = kmeans(transpose(kWTall),nClustWT,'replicates',100);
    for i = 1:nClustWT
      p(i) = loglog(centroidsWT(i,1),centroidsWT(i,2),'*','color',col(i,:));
      hold on
      idx = find(idxWT == i);
      plot(kWTall(1,idx),kWTall(2,idx),'.','color',col(i,:));
      pLeg{i} = sprintf('Cluster %d', i);
    end
    xlabel('Inflection point')
    ylabel('Slope')
    legend(p,pLeg,'location','west');

    clear p pLeg
    fName = 'FIGS/Dan-segregation-data-parameter-clustering-WT.pdf';
    saveas(gcf,fName,'pdf')
    
    
    figure('visible',visibleFlag)
    for i = 1:numel(kWTage)
      label{i} = sprintf('P%d',kWTage(i));
    end
    
    dWT = pdist(transpose(kWTall));
    LWT = linkage(dWT);
    dendrogram(LWT,0,'labels',label);
    fName = 'FIGS/Dan-segregation-data-parameter-dendrogram-WT.pdf';
    saveas(gcf,fName,'pdf')

    figure('visible',visibleFlag)
    subplot(2,1,1)
    
    p = [];
    pLeg = {};
    
    for i = 1:numel(age)
      if(~isempty(kB2KOJK{i}))
        p(end+1) = loglog(kB2KOJK{i}(1,:),kB2KOJK{i}(2,:),'.','color',col(i,:));      
        hold on
        pLeg{end+1} = sprintf('P%d',age(i));

        nB2KO(i) = size(kB2KOJK{i},2);
        kB2KOall = [kB2KOall,kB2KOJK{i}];
      else
        nB2KO(i) = 0;
      end
      
    end
    xlabel('Inflection point')
    ylabel('Slope')
    legend(p,pLeg,'location','southeast');
    title('B2KO')
    
    % We had one serious outlier
    % ax = [0 1e4 0 100]
    % axis(ax);
    
    clear p pLeg
    
    subplot(2,1,2)
    nClustB2KO = 3;
    
    [idxB2KO, centroidsB2KO] = kmeans(transpose(kB2KOall),nClustB2KO,'replicates',100);
    for i = 1:nClustB2KO
      p(i) = loglog(centroidsB2KO(i,1),centroidsB2KO(i,2),'*','color',col(i,:));
      hold on
      idx = find(idxB2KO == i);
      plot(kB2KOall(1,idx),kB2KOall(2,idx),'.','color',col(i,:));
      pLeg{i} = sprintf('Cluster %d', i);    
    end
    
    xlabel('Inflection point')
    ylabel('Slope')
    legend(p,pLeg,'location','southeast');

    % We had one serious outlier
    % axis(ax);

    fName = 'FIGS/Dan-segregation-data-parameter-clustering-B2KO.pdf';
    saveas(gcf,fName,'pdf')
    
    
    %disp(['Setting axis manually for B2KO, since we got one serious ' ...
    %      'outlier'])
    
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function plotParamSpaceALT()
try
    
    col = [228, 26, 28;
           55, 126, 184;
           77, 175, 74;
           152, 78, 163;
           255, 127, 0;
           255, 255, 51;
           166, 86, 40;
           247, 129, 191;
           153, 153, 153]/255;
    
    col = [255, 255, 217;
           237, 248, 217; 
           199, 233, 180; 
           127, 205, 187; 
           65, 182, 196; 
           29, 145, 192; 
           34, 94, 168; 
           37, 52, 148; 
           8, 29, 88]/255;
    
    nWT = zeros(numel(age),1);
    nB2KO = zeros(numel(age),1);

    yCross = [0.75 0.9];
    
    p = []; pLeg = {};
    figure('visible',visibleFlag)
    subplot(2,1,1)
    hold on

    xCrossALL = [];
    for i = 1:numel(age)
      xCross = [];
      for j = 1:numel(yCross)
        xCross(:,j) = inverseLogistic(transpose(kWTJK{i}),yCross(j));
      end
      
      try
      xCrossALL = [xCrossALL; xCross, transpose(kWTJK{i}(2,:))];
      catch e
        getReport(e)
        keyboard
      end
      
      p(i) = plot3(xCross(:,1),xCross(:,2),transpose(kWTJK{i}(2,:)),'.','color',col(i,:));
      pLeg{i} = sprintf('P%d',age(i));

    end

    xlabel(sprintf('Distance for %.2f segregation', yCross(1)))
    ylabel(sprintf('Distance for %.2f segregation', yCross(2)))    
    title('Wildtype')
    
    subplot(2,1,2)
    hold on
    p = []; pLeg = {};
    nClustWT = 5;

    [idxWT, centroidsWT] = kmeans(xCrossALL,nClustWT,'replicates',100);
    for i = 1:nClustWT
      p(i) = plot3(centroidsWT(i,1),centroidsWT(i,2),centroidsWT(i,3),'*','color',col(i,:));
      hold on
      idx = find(idxWT == i);
      plot3(xCrossALL(idx,1),xCrossALL(idx,2),xCrossALL(idx,3),'.','color',col(i,:));
      pLeg{i} = sprintf('Cluster %d', i);
    end

    xlabel(sprintf('Distance for %.2f segregation', yCross(1)))
    ylabel(sprintf('Distance for %.2f segregation', yCross(2)))    
    zlabel('Slope')
    
catch e
  getReport(e)
  keyboard
end
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
end
