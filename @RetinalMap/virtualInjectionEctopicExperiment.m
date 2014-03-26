% Retina -> SC

function [hasEctopic,fig] = virtualInjectionEctopicExperiment(obj,RGCnt, RGCdv, ...
                                                              injectionRadius, plotFlag)

  fig = []; 
  
  if(~exist('plotFlag'))
    plotFlag = 1;
  end
  
  if(~exist('injectionRadius'))
    % Check if injection site is ok
    retinaWidth = max(obj.RGCnt)-min(obj.RGCnt);
    injectionRadius = retinaWidth/2 * sqrt(0.01); % Retinal injection size
  
    % Reber (email 10 oct 2012) estimate 1-3% of retinal area
    % marked. If we assume that the retina is a disk, 
    % 0.01*pi*R^2 = pi*r^2 (r = inj radius, R = retinal radius)
    % r = R*sqrt(0.01), ie r = 0.05

  end

  debugFlag = false; % Affects colouring of clusters
  
  minClusterSize = 5; 
  distThresh = 1.5;  % Used to decide if one or two clusters
                     %
                     % !!! DW value 1.1 does not seem to work very well with Kolakov classic
  
  % 1. Find the the RGC neurons within the injection site

  r = sqrt((obj.RGCnt-RGCnt).^2 + (obj.RGCdv-RGCdv).^2);
  RGCidx = find(r <= injectionRadius);

  % Simulations with few neurons run into this problem
  if(isempty(RGCidx))
    fprintf('No neurons close to injection site, moving injection.\n')
    [~,idx] = min(r);
    RGCnt = obj.RGCnt(idx);
    RGCdv = obj.RGCdv(idx);

    r = sqrt((obj.RGCnt-RGCnt).^2 + (obj.RGCdv-RGCdv).^2);
    RGCidx = find(r <= injectionRadius);
  end
  
  
  % 2. Find the SC neurons projected to
  obj.convertConnectionTables('pre2post');
  SCidx = [];
  
  for i = 1:numel(RGCidx)
    n = obj.numPostsynapticConnections(RGCidx(i));
    SCidx = [SCidx; obj.postsynapticConnections(1:n,RGCidx(i))];
  end
  
  SCidx = unique(SCidx);
  
  
  % 3. Determine if there is an ectopic or not, using k-means clustering

  X = [obj.SCap(SCidx), obj.SCml(SCidx)];

  if(0)
    disp('Warning using GMM - less robust than K-means')
    
    % GMM turned out to be more noise sensitive when I tried it on
    % surrogate data, quite often it would just pick up the
    % spurious points I added with the second gaussian, rather than
    % separate the main cluster into two parts.
    
    % Use gaussian mixture model for ectopic detection
    options = statset('display','final');
    gm = gmdistribution.fit(X,2,'SharedCov',false,'options',options);
  
    idx = cluster(gm,X);
    SCidx1 = SCidx(find(idx == 1));
    SCidx2 = SCidx(find(idx == 2));
    
    distCent = sqrt(sum((gm.mu(1,:)-gm.mu(2,:)).^2));
    
    % Sigma = co-variances
    s1 = sqrt(gm.Sigma(1,1,1) + gm.Sigma(2,2,1));
    s2 = sqrt(gm.Sigma(1,1,2) + gm.Sigma(2,2,2));
    
    plotInjection(RGCidx,SCidx1,SCidx2);
    hold on
    if(debugFlag)
      plot(gm.mu(1,1), gm.mu(1,2),'bv', ...
           gm.mu(2,1), gm.mu(2,2),'bv')
    end
      
    if(distCent > (s1 + s2) ...
        & numel(SCidx1) >= minClusterSize & numel(SCidx2) >= minClusterSize)
      hasEctopic = true;
      title('GMM: Ectopic projection')      
    else
      hasEctopic = false;
      title('GMM: Single termination zone')      
    end
    
  else
    % Alternatively, use the old k-means (assumes same std for both regions)
    [idx,C] = kmeans(X,2,'emptyaction','singleton','replicates',5);
  
    SCidx1 = SCidx(find(idx == 1));
    SCidx2 = SCidx(find(idx == 2));
  
    % Require the clusters to contain at least 5 SC
    % The ectopic is valid if the difference of the two means is at
    % least the 1.1 * sum of the two standard deviations of the clusters.
  
    SCcent1 = C(1,:);
    SCcent2 = C(2,:);
    
    centToCentDist = sqrt(sum((SCcent1 - SCcent2).^2));
    
    % Distance to the centroid for each point in respective cluster
    centDist1 = sqrt((SCcent1(1)-obj.SCap(SCidx1)).^2 ...
                     + (SCcent1(2)-obj.SCml(SCidx1)).^2);
    centDist2 = sqrt((SCcent2(1)-obj.SCap(SCidx2)).^2 ...
                     + (SCcent2(2)-obj.SCml(SCidx2)).^2);
   
    if(plotFlag)
      plotInjection(RGCidx,SCidx1,SCidx2);
      hold on
      if(debugFlag)
        plot(SCcent1(1),SCcent1(2),'bv', ...
             SCcent2(1),SCcent2(2),'bv')
        
        fprintf('centToCent / std = %.4f\n', centToCentDist/(mean(centDist1) + mean(centDist2)))
      end
  
      if(centToCentDist > distThresh*(mean(centDist1) + mean(centDist2)) ...
         & numel(SCidx1) >= minClusterSize & numel(SCidx2) >= minClusterSize)
        hasEctopic = true;
        title('k-means: Ectopic projection')
      else
        hasEctopic = false;
        title('k-means: Single termination zone')
      end

      fName = sprintf('%s/%s-virtual-injection-retina-NT-%.2f-DV-%.2f-radius-%.3f.pdf',...
                      obj.figurePath, obj.simName, RGCnt(1), RGCdv(1),injectionRadius);
      
      switch(obj.plotFigures)
        case {0,1}
          saveas(fig,fName,'pdf');
        case 2
          obj.addParametersToFigure(strrep(fName,'.pdf','.eps'));
      end
      
    end
      
    % We get a problem if the two ectopic regions are of different
    % size, kmeans assume both clusters have same size
    
    
  end
    
  
  
  if(0)
    
    % Lets see if hierarchical clustering is better
    d = pdist(X,'euclidean');
    clustTree = linkage(d,'average');
    cophenet(clustTree,d);
    figure
    [h,nodes] = dendrogram(clustTree,0);

  end
  
  % keyboard
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function plotInjection(RGCidx, SCidx1, SCidx2)

    RGCnoIdx = setdiff(1:obj.nRGC,RGCidx);
    SCnoIdx = setdiff(1:obj.nSC,union(SCidx1,SCidx2));
  
    if(obj.plotFigures == 0)
      disp('virtualExperimentEctopic: plotFigures = 0, hiding figures.')
      visFlag = 'off';
    else
      visFlag = 'on';
    end
    
    fig = figure('visible',visFlag);

    subplot(1,2,1)
    hold on
    plot(obj.RGCnt(RGCnoIdx),obj.RGCdv(RGCnoIdx), ...
         '.', 'markersize',6, 'color', [1 1 1]*0.7);
    plot(obj.RGCnt(RGCidx),obj.RGCdv(RGCidx),'r.','markersize',12)
    
    % xlabel('Temporal - Nasal')
    % ylabel('Ventral - Dorsal')
    
    set(gca, 'xdir','reverse','ydir','reverse');
    axis equal
    mx = max(max(obj.RGCnt),1);
    my = max(max(obj.RGCdv),1);
    axis([0 mx 0 my])
        
    set(gca,'xtick',[0.1 0.9]*mx,'xticklabel',{'N','T'});
    set(gca,'ytick',[0.1 0.9]*my,'yticklabel',{'D','V'});
    set(gca,'ticklength', [0 0]);
    set(gca,'fontsize',30)
    box off
    
    subplot(1,2,2)
    hold on
    if(debugFlag)
      plot(obj.SCap(SCnoIdx),obj.SCml(SCnoIdx),'ko', ...
           obj.SCap(SCidx1),obj.SCml(SCidx1),'r*', ...
           obj.SCap(SCidx2),obj.SCml(SCidx2),'y*')    
    else
      plot(obj.SCap(SCnoIdx),obj.SCml(SCnoIdx), ...
           '.', 'markersize',6, 'color', [1 1 1]*0.7);
      plot(obj.SCap(SCidx1),obj.SCml(SCidx1),'r.', ...
           obj.SCap(SCidx2),obj.SCml(SCidx2),'r.', ...
           'markersize', 12)    
      
    end
    %xlabel('Anterior - Posterior')
    %ylabel('Medial - Lateral')
    axis equal
    mx = max(max(obj.SCap),1);
    my = max(max(obj.SCml),1);
    axis([0 mx 0 my])

    set(gca,'xtick',[0.1 0.9]*mx,'xticklabel',{'A','P'});
    set(gca,'ytick',[0.1 0.9]*my,'yticklabel',{'M','L'});
    set(gca,'ticklength', [0 0]);
    set(gca,'fontsize',30)
    box off
    
  end
  
end