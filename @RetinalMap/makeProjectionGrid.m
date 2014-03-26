% This function places a grid on the retina, and finds the
% corresponding projections in the SC.

function [grid,fig] = makeProjectionGrid(obj, nPoints, radius, debugFigs)

  fig = [];
  
  if(obj.plotFigures == 0)
    disp('superposedProjection: plotFigures = 0, hiding figures.')
    visFlag = 'off';
  else
    visFlag = 'on';
  end
  
  if(~exist('debugFigs') | isempty(debugFigs))
    debugFigs = 1;  
  end
  
  ectopicThreshold = 2;
  
  % profile on
  
  % Currently only works with eye disks
  assert(strcmp(obj.eyeType,'disk'));

  if(~exist('nPoints') | isempty(nPoints))
    nPoints = 100;
  end
  
  if(~exist('radius') | isempty(radius))
    radius = 0.05;
  end
  

  [RGCcentAP,RGCcentML] = obj.RGCprojectionCentroids();

  % 1. Randomly place a set of grid points on the retina
  % This is done by placing an abundance of putative points, and
  % then removing those that are closest together

  % Avoid grid points too close to the edge (hence padded version),
  % those grid points would have a skewed location of the centroid
  % leading to a higher density of points at edge
  % [gridX,gridY] = makeGrid(@obj.filterRetinalPositionsDisk);
  [gridX,gridY] = makeGrid(@filterRetinalPositionsDiskPadded);  

  
  % 2 Find projections for all the grid points (and their surroundings)
  % proj is n x 4 matrix, first two colums are X,Y for the major
  % projection. Last two columns are NaN, or minor ectopic.

  switch(obj.phenotype)
    case {'Isl2homozygous','Isl2heterozygous'}

      % Analayse Isl2+ and Isl2- RGCs separately
      idxP = obj.Isl2PositiveRGC;
      idxM = setdiff(1:obj.nRGC,obj.Isl2PositiveRGC);
      
      projP = findProjections(gridX,gridY,obj.RGCnt(idxP),obj.RGCdv(idxP),RGCcentAP(idxP),RGCcentML(idxP));
      projM = findProjections(gridX,gridY,obj.RGCnt(idxM),obj.RGCdv(idxM),RGCcentAP(idxM),RGCcentML(idxM));
    
    otherwise
      % Default mode
      proj = findProjections(gridX,gridY,obj.RGCnt,obj.RGCdv,RGCcentAP,RGCcentML);
  end
  
  % 3. Calculate crossovers

  % Triangulate but remove the edges that have angle smaller than angThresh
  angThresh = 20/180*pi;
  tri = makeTriangulation(angThresh);
  
  switch(obj.phenotype)
    case {'Isl2homozygous','Isl2heterozygous'}

      % In case we got Isl2 phenotype, and want to have two
      % separate maps for Isl2+ and Isl2-
      
      % 4. Calculate largest ordered submap
      [subEdgesP,numCrossingsP] = ...
          findLargestUncrossedSubGraph(tri,projP(:,1),projP(:,2),gridX,gridY);

      [subEdgesM,numCrossingsM] = ...
          findLargestUncrossedSubGraph(tri,projM(:,1),projM(:,2),gridX,gridY);
      
      
      % 5. Export the grid and the relevant data
  
      grid.RGCnt = gridX;
      grid.RGCdv = gridY;

      grid.SCap = [projP(:,1), projM(:,1)];
      grid.SCml = [projP(:,2), projM(:,2)];
      
      grid.tri = tri;
      
      grid.SCnumCrossings = [numCrossingsP, numCrossingsM];      
      
      grid.SCuncrossedGraphEdges = {subEdgesP,subEdgesM};
      
      grid.SCnodesLeft = [numel(unique(subEdgesP)), numel(unique(subEdgesM))];      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % Plot figures
      
      fig(1) = figure('visible',visFlag);
      plotProjectionIsl2(subEdgesP,subEdgesM, tri, gridX, gridY, ...
                         projP(:,1), projP(:,2), projM(:,1), projM(:,2));
    
      title(sprintf(['Largest uncrossed submap (%d, %d nodes)\n' ...
                     '(%.0f ; %.0f original crossings)'], ...
                    numel(unique(subEdgesP)), numel(unique(subEdgesM)), ...
                    numCrossingsP, numCrossingsM))
    
      figName = sprintf('%s/%s-largest-uncrossed-submap.pdf', ...
                        obj.figurePath, obj.simName);

      switch(obj.plotFigures)
        case {0,1}
          saveas(fig,figName,'pdf');
        case 2
          obj.addParametersToFigure(strrep(figName,'.pdf','.eps'));
      end    
      
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    otherwise
      % Default method
      
      % 4. Calculate largest ordered submap
      [subEdges,numCrossings] = ...
          findLargestUncrossedSubGraph(tri,proj(:,1),proj(:,2),gridX,gridY);
  
      % 5. Export the grid and the relevant data
  
      grid.RGCnt = gridX;
      grid.RGCdv = gridY;
      grid.SCap = proj(:,1);
      grid.SCml = proj(:,2);
      
      grid.tri = tri;
      grid.SCnumCrossings = numCrossings;
      grid.SCuncrossedGraphEdges = subEdges;
      grid.SCnodesLeft = numel(unique(subEdges));
  
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % Plot figure
      
      fig(1) = figure('visible',visFlag);
      plotProjection(subEdges, tri, gridX, gridY, proj(:,1), proj(:,2));
    
      title(sprintf(['Largest uncrossed submap (%d nodes)\n' ...
                     '(%.0f original crossings)'], ...
                    numel(unique(subEdges)), ...
                    numCrossings))
    
      figName = sprintf('%s/%s-largest-uncrossed-submap.pdf', ...
                        obj.figurePath, obj.simName);

      switch(obj.plotFigures)
        case {0,1}
          saveas(fig,figName,'pdf');
        case 2
          obj.addParametersToFigure(strrep(figName,'.pdf','.eps'));
      end    
    
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      if(0)
        [grid.patchEdges,grid.patchID] = locatePatches(gridX,gridY,proj,tri);
      else
        disp('Not doing patch size - code temporarilly excluded')    
      end
  end
  % profview
  
  % keyboard
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  function [NT,DV,NTpad,DVpad,mask] = filterRetinalPositionsDiskPadded(NT,DV)

    NTc = 0.5; DVc = 0.5; r = 0.5;

    % radius is the padding radius, same as the uptake region
    % radius for the grid points. We want to make sure that points
    % close to the edge are excluded. Otherwise the density of grid
    % points gets skewed, so there are more points at periphery.
    % (Comment from DS)
    
    okPos = sqrt((NT-NTc).^2+(DV-DVc).^2) < r - radius;

    idxPad = find(~okPos);
    NTpad = NT(idxPad);
    DVpad = DV(idxPad);

    idx = find(okPos);
    NT = NT(idx);
    DV = DV(idx);
  
    mask = okPos;
    
  end
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [gridX,gridY] = makeGrid(constrainFunction)
    
    % We want a fairly regular grid, so set packing factor to 0.54
    [gridX,gridY] = obj.dMin(constrainFunction,nPoints,[],0.54);
    
    gridDist = sqrt((kron(gridX,ones(1,length(gridX))) ...
                     - kron(transpose(gridX),ones(length(gridX),1))).^2 ...
                    + (kron(gridY,ones(1,length(gridY))) ...
                       - kron(transpose(gridY),ones(length(gridY),1))).^2);
     
    gridDist = gridDist + diag(diag(gridDist*inf));

    
    % Make sure that the distance between all points is at least
    % 2*radius, otherwise warn user.
    
    minSeparation = min(gridDist(:));
    if(minSeparation < 2 * radius)
      fprintf('makeProjectionGrid: minimal grid distance %f, radius %f\n',...
              minSeparation, radius)
      disp('Warning: You might want to decrease radius, or number of grid points')
    end
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Find projections for the grid points

  function [gridX,gridY] = makeGridOLD(constrainFunction)

    if(0)
      H = haltonset(2,'Skip',1e3,'Leap',1e2);
      H = scramble(H,'RR2');
      
      X0 = H(1:2*nPoints,:);
      
    else
      X0 = rand(15*nPoints,2);
    end
      
    [gridX,gridY,~,~] = constrainFunction(X0(:,1),X0(:,2));
      
    gridDist = sqrt((kron(gridX,ones(1,length(gridX))) ...
                     - kron(transpose(gridX),ones(length(gridX),1))).^2 ...
                    + (kron(gridY,ones(1,length(gridY))) ...
                       - kron(transpose(gridY),ones(length(gridY),1))).^2);
  
    gridDist = gridDist + diag(diag(gridDist*inf));
    
    [~,idx] = sort(gridDist(:));
    [idxR,idxC] = ind2sub(size(gridDist),idx);

    % Remove the grid points that are closest together until
    % nPoints remain   
  
    idxKeep = 1:numel(gridX);
    nKeep = numel(gridX);
    i = 1;
    
    while(nKeep > nPoints)
      % Find next closest pair in list, if both still remain, remove one
      if(all(ismember([idxR(i) idxC(i)],idxKeep)))
        if(rand(1) < 0.5)
          idxKeep(idxR(i)) = NaN;
        else
          idxKeep(idxC(i)) = NaN;
        end        
      
        nKeep = nKeep - 1;
      end
    
      i = i + 1;
            
    end
  
    idxKeep = find(~isnan(idxKeep));
 
    % figure
    % plot(X,Y,'o','color',[1 1 1]*0.8);
    
    gridX = gridX(idxKeep);
    gridY = gridY(idxKeep);
    
    % hold on
    % plot(X,Y,'r*')
    % hold off
    
    % Make sure that the distance between all points is at least
    % 2*radius, otherwise warn user.
    
    minSeparation = min(min(gridDist(idxKeep,idxKeep)));
    if(minSeparation < 2 * radius)
      fprintf('makeProjectionGrid: minimal grid distance %f, radius %f\n',...
              minSeparation, radius)
      disp('Warning: You might want to decrease radius, or number of grid points')
    end
          
  end
        

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % gridX,gridY = grid points
  % X,Y = coords for structure A
  % projCentX, projCentY = coords for projections to structure B
  
  function proj = findProjections(gridX,gridY,X,Y,projCentX,projCentY)
  
    minClusterSize = 5;
    
    proj = zeros(numel(gridX),4)*NaN;
    
    % 1. Find all the neurons close to selected grid point
    allProj = {};

    for i = 1:numel(gridX)
      
      d = sqrt((gridX(i)-X).^2 + (gridY(i)-Y).^2);
      idx = find(d < radius);

      if(isempty(idx))
        % Do nothing, the projections are initialised to NaN
      elseif(numel(idx) == 1)
        % Just one point
        proj(i,1:2) = [projCentX(idx), projCentY(idx)];
      else
        % More than one RGC neuron, are they clustered together?
        [group,C] = kmeans([projCentX(idx), projCentY(idx)], ...
                                2,'emptyaction','singleton','replicates',5);

        try
          idx1 = idx(find(group == 1));
          idx2 = idx(find(group == 2));
        catch e
          getReport(e)
          keyboard
        end
      
        cent1 = C(1,:);
        cent2 = C(2,:);
      
        centToCentDist = sqrt(sum((cent1 - cent2).^2));
      
        d1 = sum(sqrt((projCentX(idx1)-cent1(1)).^2 ...
                      + (projCentY(idx1)-cent1(2)).^2)) ...
             / numel(idx1);
        d2 = sum(sqrt((projCentX(idx2)-cent2(1)).^2 ...
                      + (projCentY(idx2)-cent2(2)).^2)) ...
             / numel(idx2);

        if(centToCentDist > ectopicThreshold*(d1+d2) ...
           & numel(idx1) >= minClusterSize ...
           & numel(idx2) >= minClusterSize)

          fprintf('Ectopics detected for node %.0f (C = %f)\n', ...
                  i, centToCentDist/(d1+d2))          
                    
          % We got ectopics, make the largest one the primary ectopic
          if(numel(idx1) >= numel(idx2))
            proj(i,:) = [cent1, cent2];
          else
            proj(i,:) = [cent2, cent1];
          end
          
          if(debugFigs > 1)
            figure('visible',visFlag);

            plot(projCentX(idx1),projCentY(idx1),'ko', ...
                 projCentX(idx2),projCentY(idx2),'bo')
            hold on
            p = plot(proj(i,1),proj(i,2),'r*', ...
                          proj(i,3),proj(i,4),'y*', 'markersize',12);
            hold off
            legend(p,'Major ectopic center', 'Minor ectopic center')
            axis equal
            axis([0 1 0 1])
            box off
            title(sprintf('Cent to Cent / (d1 + d2) = %f', ...
                          centToCentDist/(d1+d2)))
          end
          
        else
          % Just one projection
          proj(i,1:2) = [mean(projCentX(idx)), mean(projCentY(idx))];
          
        end
      end
    end
      
  end
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function tri = makeTriangulation(angThresh)

    tri = delaunay(gridX,gridY);

    gx = gridX(tri(:,[1 2 3 1]));
    gy = gridY(tri(:,[1 2 3 1]));   
  
    % Calculate all the edge distances
    d = sqrt(diff(gx,[],2).^2 + diff(gy,[],2).^2);

    % Get the corresponding angles
    a2 = acos((d(:,1).^2+d(:,2).^2-d(:,3).^2)./(2*d(:,1).*d(:,2)));
    a3 = acos((d(:,2).^2+d(:,3).^2-d(:,1).^2)./(2*d(:,2).*d(:,3)));
    a1 = acos((d(:,3).^2+d(:,1).^2-d(:,2).^2)./(2*d(:,3).*d(:,1)));  

    % Prune away bad angles
    goodIdx = find(a1 >= angThresh & a2 >= angThresh & a3 >= angThresh);
    tri = tri(goodIdx,:);
    
    % figure, triplot(tri,gridX,gridY);
    
  end 
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function newCount = updateIntersectionCount(oldCount,newEdges, ...
                                              removedEdges,X,Y)
    newCount = oldCount;
    
    % Loop through removed edges
    for i = 1:size(removedEdges,1)
      
      % Check for collision with remaining edges
      for j = 1:size(newEdges,1)
      
        % Skip line pairs that share one vertex                  
        if(removedEdges(i,1) ~= newEdges(j,1) & ...
           removedEdges(i,1) ~= newEdges(j,2) & ...
           removedEdges(i,2) ~= newEdges(j,1) & ...
           removedEdges(i,2) ~= newEdges(j,2))

          % Subtract the crossings from the old count to get new count
          if(doesCross(X(removedEdges(i,:)),Y(removedEdges(i,:)), ...
                       X(newEdges(j,:)),Y(newEdges(j,:)))) 
          
            nodeIdx = [removedEdges(i,:), newEdges(j,:)];
            newCount(nodeIdx) = newCount(nodeIdx) - 1;
            
          end
        end
      end
    end
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
  function [doesIntersect,P] = doesIntersect(X1,Y1,X2,Y2)
    
    % p + t*r = q + u*s
    % trick is to do cross with s
    % http://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
    % REF: "Intersection of two lines in three-space" by Ronald
    % Graham, published in Graphics Gems, page 304
      
    p = [X1(1), Y1(1), 0];
    q = [X2(1), Y2(1), 0];
    r = [X1(2)-X1(1), Y1(2)-Y1(1), 0];
    s = [X2(2)-X2(1), Y2(2)-Y2(1), 0];
    
    t = myCross(q-p,s) / myCross(r,s);
    u = myCross(q-p,r) / myCross(r,s);
    
    doesIntersect = false;
    P = [];
    
    if(myCross(r,s) == 0)
      % Lines are parallell
      
      if(myCross(q-p,r) == 0)
        % Colinear, infinite number of intersections
        doesIntersect = true;
        P = t*r;        
      end

    % Is the intersection on the line segments
    elseif(0 <= t & t <= 1 & 0 <= u & u <= 1)
      doesIntersect = true;
      P = t*r;              
    end
      
    if(0)
      % Debug plot
      figure, plot(X1,Y1,'k-',X2,Y2,'b-')
      if(doesIntersect)
        hold on
        plot(P(1),P(2),'r*')
        hold off
      end
    end
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function crossFlag = doesCross(X1,Y1,X2,Y2)
    
   % Function from Adrianna, who took it from David Sterratt,
    % who took parts of it from a webpage...

    x1x2 = X1(1) - X1(2);
    x3x4 = X2(1) - X2(2);
    y1y2 = Y1(1) - Y1(2);
    y3y4 = Y2(1) - Y2(2);
    
    dx12y12 = det([X1(1), Y1(1); X1(2), Y1(2)]);
    dx34y34 = det([X2(1), Y2(1); X2(2), Y2(2)]);
    D = det([x1x2, y1y2; x3x4,y3y4]);
    
    if(D~= 0)
      intersection_x = det([dx12y12, x1x2; dx34y34, x3x4])/D;
      intersection_y = det([dx12y12, y1y2; dx34y34, y3y4])/D;
    else
      crossFlag = false;
      return
    end
    
    % Determine if the intersection lies on the line
    on_L1 = intersection_x >= min(X1(1),X1(2)) ...
            & intersection_x <= max(X1(1),X1(2)) ...
            & intersection_y >= min(Y1(1), Y1(2)) ...
            & intersection_y <= max(Y1(1),Y1(2));
    on_L2 = intersection_x >= min(X2(1),X2(2)) ...
            & intersection_x <= max(X2(1),X2(2)) ...
            & intersection_y >= min(Y2(1), Y2(2)) ...
            & intersection_y <= max(Y2(1),Y2(2));

    if(on_L1 == 1 & on_L2 == 1)
      crossFlag = true;
    else
      crossFlag = false;
    end
    
  end
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % For each node, return the number of intersections it is
  % involved in
  
  function count = numIntersectionsPerNode(edges,X,Y)
    
    count = zeros(size(X));
    
    
    % Loop through all edges
    for i = 1:size(edges,1)
      
      % Check for collisions with all the other edges
      for j = (i+1):size(edges,1)

        % Skip line pairs that share one vertex                  
        %if(~any(intersect(edges(i,:),edges(j,:))))
        if(edges(i,1) ~= edges(j,1) & ...
           edges(i,1) ~= edges(j,2) & ...
           edges(i,2) ~= edges(j,1) & ...
           edges(i,2) ~= edges(j,2))
          
          %if(doesIntersect(X(edges(i,:)),Y(edges(i,:)), ...
          %                 X(edges(j,:)),Y(edges(j,:))))
    
          if(doesCross(X(edges(i,:)),Y(edges(i,:)), ...
                       X(edges(j,:)),Y(edges(j,:))))
            
            % Lines intersect, increment counter for the
            % corresponding nodes
            nodeIdx = [edges(i,:), edges(j,:)];
            count(nodeIdx) = count(nodeIdx) + 1;
            
          end
          
        end
        
      end
    
    end
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % edges are unique edges
  
  function [largestSubSize, subGraphID] = findSubgraphs(edges,X,Y)
  
    subGraphID = zeros(size(X));
    nextID = 1;
  
    for i = 1:size(edges,1)
      n1 = edges(i,1);
      n2 = edges(i,2);
      
      if(subGraphID(n1) == 0 & subGraphID(n2) == 0)
        % Seed a new subgraph
        subGraphID(n1) = nextID;
        subGraphID(n2) = nextID;
        nextID = nextID + 1;
             
      elseif(subGraphID(n1) == subGraphID(n2))
        % Already connected do nothing
        
      elseif(subGraphID(n1) == 0)
        % Add to existing subgraph (we know both cant be zero)
        subGraphID(n1) = subGraphID(n2);
        
      elseif(subGraphID(n2) == 0)
        % Add to existing subgraph
        subGraphID(n2) = subGraphID(n1);
        
      else
        % Both are non-zero, merge the two subgraphs
        subGraphID(find(subGraphID == subGraphID(n2))) = subGraphID(n1);
        
      end
        
    end
    
    % Find the size of the largest connected component
    % 0 marks deleted nodes, not included
    usedID = setdiff(unique(subGraphID),0);
    
    subGraphSize = zeros(size(usedID));
    
    for i = 1:numel(usedID)
      subGraphSize(i) = nnz(subGraphID == usedID(i));
    end
      
    largestSubSize = max(subGraphSize);
    
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function numDisconnects = countDisconnects(edges,X,Y)
    
    [largestOriginalComponent, ~] = findSubgraphs(edges,X,Y);
    
    % For each node, see how many nodes would be disconnected if we
    % remove that node
    
    numDisconnects = zeros(size(X));
    connectedNodes = unique(edges(:));
      
    for i = 1:numel(connectedNodes)
      
      [removeIdx,foo] = find(edges == connectedNodes(i));
      edgeIdx = setdiff(1:size(edges,1),removeIdx);
      
      [largestComp, subGraphID] = findSubgraphs(edges(edgeIdx,:),X,Y);
      numDisconnects(connectedNodes(i)) = ...
          largestOriginalComponent - largestComp;
            
    end
          
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function [edges,numCrossings] = ...
        findLargestUncrossedSubGraph(tri,X,Y,gridX,gridY)
  
    % gridX and gridY are just for plotting purposes
    
    nodeIdx = unique(tri(:));
    done = false;

    % Extract all the unique edges
    edges = [tri(:,1), tri(:,2); 
             tri(:,2), tri(:,3); 
             tri(:,3), tri(:,1)];
    edges = unique(sort(edges,2),'rows');

    count = numIntersectionsPerNode(edges,X,Y);    
    
    % Each crossing is counted at both nodes of the two edges that
    % cross, ie 4 times
    numCrossings = sum(count/4); 
    
    if(debugFigs > 1)
      plotEdges(edges,X,Y)
      hold on
      for i = 1:numel(X)
        text(X(i),Y(i),sprintf('%.0f',count(i)))
      end
      hold off
      box off
      title(sprintf('Total crossings %.0f', numCrossings))
    end
    
    while(~done)
      
      if(0)
        % Test verification
        testCount = numIntersectionsPerNode(edges,X,Y);
        if(any(testCount ~= count))
          disp('I fucked up!')
          keyboard
        end
      end
      
      numDisconnects = countDisconnects(edges,X,Y);

      % Debug plot
      if(0)
        figure
        subplot(2,2,1)
        plotEdges(edges,X,Y)
        hold on
        for i = 1:numel(numDisconnects)
          text(X(i),Y(i),sprintf('%.0f',numDisconnects(i)));
        end
        hold off
        title('Num disconnects')
        subplot(2,2,2)
        plotEdges(edges,X,Y)
        hold on
        for i = 1:numel(count)
          text(X(i),Y(i),sprintf('%.0f',count(i)));
        end
        title('Num crossings')
        hold off
        subplot(2,2,3)
        plotEdges(edges,X,Y)
        hold on
        for i = 1:numel(X)
          text(X(i),Y(i),sprintf('%.0f',i),'color',[1 0 0]);
        end
        hold off
        title('Node number')   
        
        % keyboard
      end
            
      % We want to remove the node that has the most cross-overs,
      % while not disconnecting more nodes than itself.
      candidateNodes = find(numDisconnects == 1 & count > 0);
    
      if(isempty(candidateNodes))
        done = true;
      else
        % Remove the one with the most crossings of the candidates
        [~,maxIdx] = max(count(candidateNodes));

        killNodeIdx = candidateNodes(maxIdx);
        % fprintf('Removing node %d\n', killNodeIdx)

        [rowIdx,colIdx] = find(edges == killNodeIdx);
        
        removedEdges = edges(rowIdx,:);
        edges(rowIdx,:) = [];
        
        count = updateIntersectionCount(count,edges,removedEdges,X,Y);
      
        [edges, count] = pruneGridTails(edges, count, X, Y);
        
      end
      
    end

    % Remove the tail bits
    edges = pruneGridTails(edges, count, X, Y);    
        
    % keyboard
    
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Faster less general version of the cross function in matlab
  
  function c = myCross(a,b)
    
    c = [a(2)*b(3)-a(3)*b(2),a(1)*b(3)-a(3)*b(1),a(1)*b(2)-a(2)*b(1)];
    
    % assert(all(c == cross(a,b)))
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function plotEdges(edges,X,Y)
    
    plot(X(transpose(edges)),Y(transpose(edges)), ...
         '-','color',[1 1 1]*0.7);
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Remove the tails  
  
  function [prunedEdges, newCount] = pruneGridTails(edges, count,X,Y)
    
    % Find the nodes that only have one edge, remove them, then
    % repeat until all remaining nodes have two or more edges
      
    % disp('Pruning away the tails from the grid.')
      
    prunedEdges = edges;
    nRemoved = -1;
    
    newCount = count;
    
    while(nRemoved ~= 0)
      nRemoved = 0;
      
      A = accumarray(prunedEdges(:),1);
      singleNodes = find(A == 1);
      
      for i = 1:numel(singleNodes)

        % disp('Pruneing tails...')
        
        % Remove the edges that connect the single nodes
        [rowIdx,colIdx] = find(prunedEdges == singleNodes(i));
        
        removeEdges = prunedEdges(rowIdx,:);
        
        prunedEdges(rowIdx,:) = [];
        nRemoved = nRemoved + numel(rowIdx);
        
        newCount = updateIntersectionCount(newCount,prunedEdges, ...
                                           removeEdges, X, Y);
      
      end
      
    end
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function [patchEdges,patchID] = locatePatches(gridX,gridY,proj,tri)
try
    
    % Extract all the unique edges
    edges = [tri(:,1), tri(:,2); 
             tri(:,2), tri(:,3); 
             tri(:,3), tri(:,1)];
    edges = unique(sort(edges,2),'rows');

    X = proj(:,1);
    Y = proj(:,2);
    
    % Sort all edges according to how much they are stretched 
    % after projection
    rDist = sqrt((gridX(edges(:,1))-gridX(edges(:,2))).^2 ...
                 + (gridY(edges(:,1))-gridY(edges(:,2))).^2);
    
    scDist = sqrt((X(edges(:,1))-X(edges(:,2))).^2 ...
                  + (Y(edges(:,1))-Y(edges(:,2))).^2);
    
    edgeStretch = scDist./rDist;
    
    count = numIntersectionsPerNode(edges,proj(:,1),proj(:,2));

    [eStretch, sortIdx] = sort(edgeStretch,'descend');
    
    sortedEdges = edges(sortIdx,:);
    nextIdx = 1;
    
    % Repeat if there are still crossings        
    while(sum(count) > 0)
      
      % Candidate for removal, the most stretched edge
      removeEdge = sortedEdges(nextIdx,:);
      nextIdx = nextIdx + 1;
      
      if(any(count(removeEdge)) > 0)
        idx = find(edges(:,1) == removeEdge(1) ...
                   & edges(:,2) == removeEdge(2));
        edges(idx,:) = [];
            
        % Update the count
        count = updateIntersectionCount(count,edges,removeEdge,X,Y);
      else
        
        % Neither of the two nodes connected by that edge
        % participated in any crossings, keep edge. 
        
      end
        
    end
      
    % All crossings are gone, find connected patches
    patchID = zeros(numel(gridX),1);
    nextID = 1;
    
    for i = 1:size(sortedEdges,1)
      if(~any(patchID(sortedEdges(i,:))))
        % Both nodes are unconnected, create new patch
        patchID(sortedEdges(i,:)) = nextID;
        nextID = nextID + 1;
      elseif(all(patchID(sortedEdges(i,:))))
        % Both are member of a patch, merge them
        patchID(find(patchID == patchID(sortedEdges(i,1)))) = ...
            patchID(sortedEdges(i,2));
      else
        % One belongs to a patch, the other one does not, add the
        % lone node to the patch
        patchID(sortedEdges(i,:)) = max(patchID(sortedEdges(i,:)));
      end
    end
    
    % Separate sub patches that are only connected by one link
    
    % !!!!!!!
    
    patchEdges = edges;
   
    if(debugFigs)
      
      fig(2) = figure('visible',visFlag);
      plotProjection(patchEdges, tri, gridX, gridY, X, Y);
      
    end
    
catch e
  getReport(e)
  keyboard
end
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function plotProjection(edges, tri, gridX, gridY, X, Y)
    
    okIdx = unique(edges(:));
    badIdx = setdiff(1:max(tri(:)),okIdx);

    subplot(1,2,1)
    plotEdges(edges,gridX,gridY);
    hold on
    plot(gridX(okIdx),gridY(okIdx),'ko', ...
         gridX(badIdx),gridY(badIdx),'ro');
    hold off
    set(gca,'xdir','reverse','ydir','reverse');
    box off
    axis equal
    mx = max(max(obj.RGCnt),1);
    my = max(max(obj.RGCdv),1);
    axis([0 mx 0 my])
        
    set(gca,'xtick',[0.1 0.9]*mx,'xticklabel',{'N','T'});
    set(gca,'ytick',[0.1 0.9]*my,'yticklabel',{'D','V'});
    set(gca,'ticklength', [0 0]);        
    set(gca,'fontsize',30)
    
    
    title(obj.simName)
    subplot(1,2,2)
    plotEdges(edges,X,Y)
    hold on
    plot(X(okIdx),Y(okIdx),'ko', ...
         X(badIdx),Y(badIdx),'ro')
    hold off

    axis equal
    mx = max(max(obj.SCap),1);
    my = max(max(obj.SCml),1);
    axis([0 mx 0 my])
    set(gca,'xtick',[0.1 0.9]*mx,'xticklabel',{'A','P'});
    set(gca,'ytick',[0.1 0.9]*my,'yticklabel',{'M','L'});
    set(gca,'ticklength', [0 0]);
    box off
    set(gca,'fontsize',30)    
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function plotProjectionIsl2(edgesP, edgesM, tri, gridX, gridY, XP, YP, XM, YM)
    
    okIdxP = unique(edgesP(:));
    okIdxM = unique(edgesM(:));
    
    okIdx = intersect(okIdxP,okIdxM);
    partialOkIdx = union(setdiff(okIdxP,okIdxM),setdiff(okIdxM,okIdxP));
    badIdx = setdiff(1:max(tri(:)),union(okIdxP,okIdxM));
    
    badIdxP = setdiff(1:max(tri(:)),okIdxP);
    badIdxM = setdiff(1:max(tri(:)),okIdxM);

    subplot(1,2,1)
    hold on
    plotEdges(edgesP,gridX,gridY);
    plotEdges(edgesM,gridX,gridY);    

    plot(gridX(okIdx),gridY(okIdx),'ko', ...
         gridX(partialOkIdx),gridY(partialOkIdx),'ro', ...
         gridX(badIdx),gridY(badIdx),'rx');
    hold off
    box off
    set(gca,'xdir','reverse','ydir','reverse');
    
    axis equal
    
    mx = max(max(obj.RGCnt),1);
    my = max(max(obj.RGCdv),1);
    axis([0 mx 0 my])
        
    set(gca,'xtick',[0.1 0.9]*mx,'xticklabel',{'N','T'});
    set(gca,'ytick',[0.1 0.9]*my,'yticklabel',{'D','V'});
    set(gca,'ticklength', [0 0]);        
    set(gca,'fontsize',30)
 
      
    title(obj.simName)
    subplot(1,2,2)
    hold on
    plotEdges(edgesP,XP,YP)
    plotEdges(edgesM,XM,YM)
    
    plot(XP(okIdxP),YP(okIdxP),'bo', ...
         XM(okIdxM),YM(okIdxM),'ko', ...
         XP(badIdxP),YM(badIdxP),'ro', ...
         XM(badIdxM),YM(badIdxM),'ro')
    hold off

    axis equal
    mx = max(max(obj.SCap),1);
    my = max(max(obj.SCml),1);
    axis([0 mx 0 my])
    set(gca,'xtick',[0.1 0.9]*mx,'xticklabel',{'A','P'});
    set(gca,'ytick',[0.1 0.9]*my,'yticklabel',{'M','L'});
    set(gca,'ticklength', [0 0]);
    box off
    set(gca,'fontsize',30)  
    
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
end
