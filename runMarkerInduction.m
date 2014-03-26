function obj=runMarkerInduction(obj, method)
  if (~exist('method')) 
    method = 'C';
  end
  
  try
    % Mapping between quantities Willshaw (2006) and Johannes' framework and
    % 
    % RAi - obj.RGCEphA          
    % RBi - obj.RGCEphB          
    % TAj - obj.SCephrinA 
    % TBj - obj.SCephrinB 
    % Sij - obj.connectionMatrix
    % IjA - Calculated in stepMarkerInduction()
    % IjB - Calculated in stepMarkerInduction()
    
    % Use the Delaunay triangulation to find nearest Neighbours
    % figure
    tSC = delaunay(obj.SCap, obj.SCml);
    % trimesh(tSC, obj.SCap, obj.SCml, 'Color', 'r');    
    % hold on
    % Filter out triangles that have angles of 10 degrees or less
    R = [obj.SCap obj.SCml];            % Coords
    tSC = tSC(angles(tSC(:, [1 2 3]), R) > 10,:);
    tSC = tSC(angles(tSC(:, [2 3 1]), R) > 10,:);
    tSC = tSC(angles(tSC(:, [3 1 2]), R) > 10,:);
    % trimesh(tSC, obj.SCap, obj.SCml, 'Color', 'k');
    
    assert(all(size(obj.SCap) == size(obj.SCephrinA)))
    d = obj.SCDistance();

    % Neighbours matrix method of determining neighbours
    neighbours = zeros(obj.nSC, obj.nSC);
    %         set up links between  neighbours
    for i=1:size(tSC,1);
      neighbours(tSC(i,1), tSC(i,2)) = 1;
      neighbours(tSC(i,2), tSC(i,1)) = 1;
      neighbours(tSC(i,1), tSC(i,3)) = 1;
      neighbours(tSC(i,3), tSC(i,1)) = 1;
      neighbours(tSC(i,2), tSC(i,3)) = 1;
      neighbours(tSC(i,3), tSC(i,2)) = 1;
    end
    
    % For each SC cell, find indicies of all connected cells
    for i = 1:obj.nSC
      idx = find(neighbours(:, i) == 1);
      assert(length(idx) == length(unique(idx)))
      obj.neighbourSC{i} = int32(idx);
      obj.nNeighbourSC(i) = numel(idx);
      obj.info.neighbourSCdist{i} = d(idx, i);
    end
    obj.info.IA = zeros(obj.nSC, 1);
    obj.info.IB = zeros(obj.nSC, 1);

    if (obj.plotFigures == 0)
      disp('runMarkerInduction: plotFigures = 0, hiding figures!')
      visFlag = 'off';
    else
      visFlag = 'on';
    end
    fig = figure('visible', visFlag);
    
    % Find times at which to report iterations
    tN = 60000/obj.info.dt;
    t0 = obj.reportStep;
    N = 10; 
    reportIterSteps = t0*(tN/t0).^((0:(N-1))./(N-1))
    reportIterSteps = sort([reportIterSteps 40000/obj.info.dt 50000/obj.info.dt 70000/obj.info.dt ...
                       80000/obj.info.dt 90000/obj.info.dt 100000/obj.info.dt ...
                       110000/obj.info.dt 120000/obj.info.dt])
    itrep = 1;                          % Index into reportIterSteps
    
    % Now run the simulation
    cputime0 = cputime;
    while(obj.curStep < obj.nSteps)
      if (obj.plotFigures == 1)
        % Plotting (if desired)
        subplot(331)
        plot(obj.SCap, obj.info.IA, 'r.')
        title('Induced EphA')
        subplot(334)
        plot(obj.SCml, obj.info.IB, 'b.')
        title('Induced EphB')
        subplot(332)
        plot(obj.SCap, obj.SCephrinA, 'r.')
        title({['t = ' num2str(obj.curStep*obj.info.dt)];'ephrinA'})
        subplot(335)
        plot(obj.SCml, obj.SCephrinB, 'b.')
        title('ephrinB')
        subplot(333)
        plot(obj.SCap, 1 - obj.info.RGCEphAScale .* obj.SCephrinA .* obj.info.IA, 'r.')
        title('1 - RGCEphAScale T^A I^A')
        subplot(336)
        plot(obj.SCml, obj.info.IB - obj.SCephrinB, 'b.')
        title('I^B - T^B')
        
        % Find projection of SC cells on Retina
        % Find RGC that maximally excites each SC cell
        [S, I] = max(obj.connectionMatrix, [], 2);
        DT = delaunay(obj.SCap, obj.SCml);
        SCapProj = obj.RGCnt(I);
        SCmlProj = obj.RGCdv(I);
        subplot(337)
        triplot(DT, SCapProj, SCmlProj);

        % Plot colourmaps of ephrinA and ephrinB
        cm = colormap(jet);
        subplot(338)
        cols = interp1(linspace(min(obj.SCephrinA), ...
                                max(obj.SCephrinA+1e-6), ...
                                size(cm, 1)), ...
                       cm, obj.SCephrinA);
        cla
        hold on
        for i=1:obj.nSC
          plot(obj.SCap(i), obj.SCml(i), '.', 'Color', cols(i,:))
        end
        title('ephrinA')

        subplot(339)
        cols = interp1(linspace(min(obj.SCephrinB), ...
                                max(obj.SCephrinB), ...
                                size(cm, 1)), ...
                       cm, obj.SCephrinB);
        cla
        hold on
        for i=1:obj.nSC
          plot(obj.SCap(i), obj.SCml(i), '.', 'Color', cols(i,:))
        end
        title('ephrinB')
        
        drawnow;
      end
      assert(all(size(obj.SCap) == size(obj.SCephrinA)))
      iterSteps = min(obj.reportStep, obj.nSteps - obj.curStep);
      elapsedTime = cputime - cputime0;
      disp(['Simulating to step ' num2str(obj.curStep) ' of ' ...
            num2str(obj.nSteps) ' Elapsed time = ' num2str(elapsedTime) ...
            's ; ' num2str(elapsedTime/obj.curStep * (obj.nSteps - ...
                                                      obj.curStep) ...
                           ) ' to go']);
      if (method == 'matlab') 
        obj = stepMarkerInduction(obj, iterSteps);    
      else
        try
          [S TA TB IA IB] =  ...
              stepMarkerInductionC(transpose(obj.connectionMatrix), ...
                                   obj.SCephrinA, ...
                                   obj.SCephrinB, ...
                                   obj.RGCEphA, ...
                                   obj.RGCEphB, ...
                                   obj.neighbourSC, ...
                                   obj.info.alpha, ...
                                   obj.info.beta, ...
                                   obj.info.gamma, ...
                                   obj.info.kappa, ...
                                   obj.info.RGCEphAScale,...
                                   obj.info.dt, ...
                                   iterSteps);
        catch e
          getReport(e)
          disp('Is the source compiled? mex stepMarkerInduction.c')
          disp('Then try again. (or set method = matlab)')
          
          keyboard
        end
        obj.connectionMatrix = transpose(S);
        obj.SCephrinA = TA;
        obj.SCephrinB = TB;
        obj.info.IA = IA;
        obj.info.IB = IB;
        obj.curStep = obj.curStep + iterSteps;
      end
      assert(~any(isnan(obj.SCephrinA)));

      connectionMatrix0 = obj.connectionMatrix;
      % Save state
      obj.connectionMatrix(obj.connectionMatrix < 0.001) = 0;
      obj.convertConnectionTables('mat2pre');
      obj.connectionMatrix = connectionMatrix0;
      obj.saveState();
      if (obj.curStep >= reportIterSteps(itrep))
          obj.saveIter();
          itrep = itrep + 1;
      end
    end
  catch e
    getReport(e)
    keyboard
  end
end

% Compute angles between lines of a triangulation. Suppose that A is
% the vertex denoted R(T(:,1),:), B the vertex R(T(:,2) and C the
% vertex R(T(:,3)). Then the angle BAC is returned by this function.
function BAC=angles(T, R) 
  AB = R(T(:,2),:) - R(T(:,1),:);
  AC = R(T(:,3),:) - R(T(:,1),:);
  BAC = 180/pi*acos(sum(AB .* AC, 2)./...
                    (sqrt(sum(AB.^2, 2)).*sqrt(sum(AC.^2, 2))));
end