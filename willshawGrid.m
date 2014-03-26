% This function takes a RetinalMap object and returns the results
% from the Willshaw grid analysis
%

function [results,figHandles,P] = willshawGrid(obj,fakeID,projType)

  if(~exist('obj'))
    disp('Did you forget to specify the RetinalMap object?')
    results = [];
    figHandles = [];
    return
  end
  
  if(~exist('projType'))
    % projType = 'centroid';
    projType = 'maxRtoSC';    
  end
  

  if(~exist('fakeID') | isempty(fakeID))
    fakeID = 1000;
  end
   
  N = 100;
  radius = 0.07; %0.07;
  ectopicnodes = true; 

  addpath('../map-analysis/lib','-begin');  
  
  disp('Warning calling this code will change the random number generator')
  
  try
    
    


    switch(obj.phenotype)
      case {'Isl2homozygous','Isl2heterozygous'}

        
        if(0)
          % Do separate analysis for Isl2 positive and negative cells
          PP = getWillshawLattice(obj.Isl2PositiveRGC);
          PM = getWillshawLattice(setdiff(1:obj.nRGC,obj.Isl2PositiveRGC));
        
        else
          % Use the same retinal grid points
          [PP,PM] = getWillshawLatticeIsl2();
        end
        
        % Ideally these two plots would be in the same figure...
        figHandles(1) = do_plot(PP);
        print('-dpdf', [obj.figurePath, '/', obj.simName, '_lattice_orientation_P.pdf']);
        figHandles(2) = do_plot(PM);
        print('-dpdf', [obj.figurePath, '/', obj.simName, '_lattice_orientation_M.pdf']);

        f1 = do_lattice_plot(PP);
        print('-dpdf', [obj.figurePath, '/', obj.simName, '_latticeP.pdf']);
        f2 = do_lattice_plot(PM);
        print('-dpdf', [obj.figurePath, '/', obj.simName, '_latticeM.pdf']);
        
        % Also do full analysis
        % P = getWillshawLattice();                    
                        
      otherwise
        P = getWillshawLattice();
        figHandles = do_plot(P);
        print('-dpdf', [obj.figurePath, '/', obj.simName, '_lattice_orientation.pdf']);

        f1 = do_lattice_plot(P);
        print('-dpdf', [obj.figurePath, '/', obj.simName, '_lattice.pdf']);
        
    end

    
  catch e
    getReport(e);
    keyboard
    % Clear return variables to indicate that something failed!
    results = [];
    figHandles = [];

    % Return back to our original directory
    % cd(curDir);
  
    % Return to normal random generation
    obj.initializeRandomGenerator()
    rmpath('../map-analysis/lib');

    return
  end
    
  % Return back to our original directory
  % cd(curDir);
  
  % Return to normal random generation
  % rng('default')
  % rng('shuffle')

  % Extract the measures we want
  switch(obj.phenotype)
    case {'Isl2homozygous','Isl2heterozygous'}

      try
        results{1} = PP.stats;
        results{2} = PM.stats;
      
        P(1) = PM;
        P(2) = PP;
      catch e
			  getReport(e)
        keyboard
			end


    otherwise
      results = P.stats;
  end
  
  % results = P;
  
  rmpath('../map-analysis/lib');

  % Restore random seed generator (but with new seed)
  obj.initializeRandomGenerator();
  
  % keyboard
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  function flag = nonIntegerValues(mat)
    flag = any(abs(mat(:) - round(mat(:))) > 1e-9);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function P = getWillshawLattice(RGCidx)
  
    if(~exist('RGCidx') | isempty(RGCidx))
      RGCidx = 1:obj.nRGC;
      allRGC = true;
    else
      allRGC = false;
    end
    
    startDir = pwd;
    
    % cd('../WillshawAnalysis');

    % Initialise
    % Everything is done through passing of the P structure
    P = getparams(fakeID,N,radius);
    % Make sure that each circle has at least 3 points in it, but no more
    % than 10
    P.coll_min_points =  max(3, min(length(RGCidx)/N/2, 10));
    P.field_min_points = max(3, min(length(RGCidx)/N/2, 10));
    
    P.field.title = 'Retina';
    P.field.xlabel = 'Temporal-Nasal';
    P.field.ylabel = 'Ventral-Dorsal';
    
    % Temporarilly back to the RetinalMap directory
    % cd(curDir)
    
    % Should we use connection matrix or connection table?

    % Normally we should have a presynaptic connection table
    % * If the table does not exist...
    % * If the table exists, but has non-integer values,
    % and the connection matrix exists...
    % 
    % Then use the connection matrix instead

    if(isempty(obj.presynapticConnections) ...
       | (~isempty(obj.connectionMatrix) ...
          & nonIntegerValues(obj.presynapticWeight)))
      disp('Using connection matrix instead of connection table')
      useTableFlag = false;
    else
      disp('Using connection table (default)')
      useTableFlag = true;
    end
    
    switch(projType)
      case 'intrinsic'
        
        assert(allRGC); % Make sure we use all RGC
        
        [retinaNT,retinaDV, SCAP, SCML] = obj.virtualIntrinsicImagingProjections();
        
        P.full_field(:,1) = 1 - retinaNT();
        P.full_field(:,2) = 1 - retinaDV();
        
        % P.full_coll(:,1) = 1 - SCAP();   
        % P.full_coll(:,2) = 1 - SCML();
        
        P.full_coll(:,1) = SCAP();   
        P.full_coll(:,2) = SCML();
              
      case 'centroid'

        assert(allRGC); % Make sure we use all RGC        
        
        P.full_field(:,1) = 1 - obj.RGCnt();
        P.full_field(:,2) = 1 - obj.RGCdv();

        % This function should also be updated to handle useTableFlag
        [RGCcentAP,RGCcentML] = obj.RGCprojectionCentroids();
          
        %P.full_coll(:,1) = 1 - RGCcentAP;   
        %P.full_coll(:,2) = 1 - RGCcentML;

        P.full_coll(:,1) = RGCcentAP;   
        P.full_coll(:,2) = RGCcentML;
        
      case 'maxRtoSC'
        P.full_field(:,1) = 1 - obj.RGCnt(RGCidx);
        P.full_field(:,2) = 1 - obj.RGCdv(RGCidx);
        
        [RGCmaxAP,RGCmaxML] = obj.RGCprojectionMax(useTableFlag);
        
        % P.full_coll(:,1) = 1 - RGCmaxAP(RGCidx);   
        % P.full_coll(:,2) = 1 - RGCmaxML(RGCidx);

        P.full_coll(:,1) = RGCmaxAP(RGCidx);   
        P.full_coll(:,2) = RGCmaxML(RGCidx);

        
      case 'maxSCtoR'

        assert(allRGC); % Make sure we use all RGC        
        
        % P.full_coll(:,1) = 1 - obj.SCap;   
        % P.full_coll(:,2) = 1 - obj.SCml;

        P.full_coll(:,1) = obj.SCap;   
        P.full_coll(:,2) = obj.SCml;
        
        [SCmaxNT,SCmaxDV] = obj.SCprojectionMax(useTableFlag);
        
        P.full_field(:,1) = 1 - SCmaxNT;
        P.full_field(:,2) = 1 - SCmaxDV;
        
      otherwise
        fprintf('Unknown projection method %s\n', projType)
        keyboard
    end

    P.datalabel = obj.simName;
        
    try
      P = run_data(P);
    catch e
      disp('Oh no, it crashed in the lattice analysis!')
      getReport(e)
      keyboard
    end
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function [PP,PM] = getWillshawLatticeIsl2()
  
    startDir = pwd;
    
    % cd('../WillshawAnalysis');

    % Initialise
    % Everything is done through passing of the P structure
    P = getparams(fakeID, N, radius);
    
    % Temporarilly back to the RetinalMap directory
    % cd(curDir)
    
    % Should we use connection matrix or connection table?

    % Normally we should have a presynaptic connection table
    % * If the table does not exist...
    % * If the table exists, but has non-integer values,
    % and the connection matrix exists...
    % 
    % Then use the connection matrix instead

    if(isempty(obj.presynapticConnections) ...
       | (~isempty(obj.connectionMatrix) ...
          & nonIntegerValues(obj.presynapticWeight)))
      disp('Using connection matrix instead of connection table')
      useTableFlag = false;
    else
      disp('Using connection table (default)')
      useTableFlag = true;
    end
    
    switch(projType)
      case 'intrinsic'
        
        assert(allRGC); % Make sure we use all RGC
        
        [retinaNT,retinaDV, SCAP, SCML] = obj.virtualIntrinsicImagingProjections();
        
        P.full_field(:,1) = 1 - retinaNT();
        P.full_field(:,2) = 1 - retinaDV();
        
        % P.full_coll(:,1) = 1 - SCAP();   
        % P.full_coll(:,2) = 1 - SCML();
        
        P.full_coll(:,1) = SCAP();   
        P.full_coll(:,2) = SCML();
              
      case 'centroid'

        assert(allRGC); % Make sure we use all RGC        
        
        P.full_field(:,1) = 1 - obj.RGCnt;
        P.full_field(:,2) = 1 - obj.RGCdv;

        % This function should also be updated to handle useTableFlag
        [RGCcentAP,RGCcentML] = obj.RGCprojectionCentroids();
          
        %P.full_coll(:,1) = 1 - RGCcentAP;   
        %P.full_coll(:,2) = 1 - RGCcentML;

        P.full_coll(:,1) = RGCcentAP;   
        P.full_coll(:,2) = RGCcentML;
        
      case 'maxRtoSC'
        P.full_field(:,1) = 1 - obj.RGCnt;
        P.full_field(:,2) = 1 - obj.RGCdv;
        
        [RGCmaxAP,RGCmaxML] = obj.RGCprojectionMax(useTableFlag);
        
        % P.full_coll(:,1) = 1 - RGCmaxAP;   
        % P.full_coll(:,2) = 1 - RGCmaxML;

        P.full_coll(:,1) = RGCmaxAP;   
        P.full_coll(:,2) = RGCmaxML;

        
      case 'maxSCtoR'

        assert(allRGC); % Make sure we use all RGC        
        
        % P.full_coll(:,1) = 1 - obj.SCap;   
        % P.full_coll(:,2) = 1 - obj.SCml;

        P.full_coll(:,1) = obj.SCap;   
        P.full_coll(:,2) = obj.SCml;
        
        [SCmaxNT,SCmaxDV] = obj.SCprojectionMax(useTableFlag);
        
        P.full_field(:,1) = 1 - SCmaxNT;
        P.full_field(:,2) = 1 - SCmaxDV;
        
      otherwise
        fprintf('Unknown projection method %s\n', projType)
        keyboard
    end

    P.datalabel = obj.simName;

    % Change to the directory for the willshaw analysis code
    % cd('../WillshawAnalysis');
    
    
    %%% Copied from DSrun_data.m

    % For Isl2, we split this into two cases
    PM = P;
    PP = P;
    
    PP.full_coll = P.full_coll(obj.Isl2PositiveRGC,:);
    PP.full_field = P.full_field(obj.Isl2PositiveRGC,:);

    idx = setdiff(1:obj.nRGC,obj.Isl2PositiveRGC);
    PM.full_coll = P.full_coll(idx,:);
    PM.full_field = P.full_field(idx,:);

    PP = run_data(PP);
    PM = run_data(PM);
  end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function fig = do_plot(P) 
    
    if(obj.plotFigures)
      fig = figure('visible', 'on');
    else
      fig = figure('visible', 'off');
    end
    h1 = subplot2(1, 3, 1);
    h2 = subplot2(1, 3, 2);
    h3 = subplot2(1, 3, 3);
    
    plot_lattice(P, 'FTOC', h1, h2, ...
                 'Subgraph', true, ...
                 'EctOptions', 4, ...
                 'AncSize', 10, ...
                 'AxisStyle', 'none', ...
                 'Outline', 'none');
    
    % Plot SC outline
    set(gcf,'currentaxes',h2);
    idx = convhull(obj.SCap,obj.SCml);
    pBorder = plot(obj.SCap(idx),obj.SCml(idx),'color',[1 1 1]*0.7);
    set(pBorder,'z',-1*ones(numel(idx),1));
    
    % 'AxisStyle', 'box');      
    plot_angles(P, 'FTOC', h3)
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function fig = do_lattice_plot(P) 
    
    if(obj.plotFigures)
      fig = figure('visible', 'on');
    else
      fig = figure('visible', 'off');
    end
    h1 = subplot2(1, 2, 1);
    h2 = subplot2(1, 2, 2);
    
    col = [30, 142, 59; 
           0, 165, 190; 
           65, 182, 196; 
           44, 127, 184; 
           37, 52, 148; 
           254, 204, 92; 
           253, 141, 60; 
           240, 59, 32; 
           189, 0, 38] / 255;
    
    % Make centre colour look better matched to both
    col(3,:) = [0 0 0];
    
    % Permute order to match lattice
    col = col([1 2 9 8 3 7 6 4 5],:);
    
    % !!! TODO
    % 1. Increase Anchor size
    % 2. Decrease line width to 0.1 
    % 3. Change colour for x-axis?
    plot_lattice(P, 'FTOC', h1, h2, ...
                 'Subgraph', true, ...
                 'EctOptions', 4, ...
                 'LatticeLineWidth', 0.2, ...
                 'AncSize', 15, ...
                 'AncColours', col, ...
                 'AxisStyle', 'none', ...
                 'Outline', 'none');
    
    % Plot SC outline
    set(gcf,'currentaxes',h2);
    idx = convhull(obj.SCap,obj.SCml);
    pBorder = plot(obj.SCap(idx),obj.SCml(idx),'color',[1 1 1]*0.7);
    set(pBorder,'z',-1*ones(numel(idx),1));
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
end
