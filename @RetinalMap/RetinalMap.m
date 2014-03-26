classdef RetinalMap < handle

  properties

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Simulation ID and directories used
    
    simName = 'default';
    dataPath = 'SAVE';
    figurePath = 'FIGS';
    
    simID = 0;
    timeStamp = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % CSV file with gradients specified for different phenotypes
    
    gradientInfoFile = 'gradients/Eph-ephrins-both-axis-JH.csv';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    curStep = 0; % This is the current step, normally when I
                 % display iterations, I choose to do it as average
                 % number of iterations per neuron passed.
    time = 0;
    
    nSteps = 1000;
    reportStep = 100;

    nRGC = 10000;
    nSC = 10000;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Parameters for placing neurons
    
    % optimalPackingFactor = 1/6*pi*sqrt(3) = 0.9069 -- Steinhaus 1999, p 202
    dMinPackingFactor = 0.4; 
    % dMinPF = 0.4 Gives Voronoi regularity index 5.2.
    % Raven, Stagg, Nassar, Reese 2005 gives 5 (P1) to 6.5 (P10)
    dMinPackingFactor3D = 0.3;
    
    % PF = 0.28 for horisontal cells, Raven, Stagg, Nassar, Reese 2005
    % 0.5473 -- Masharu Tanemura, 1979 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    debugFlag = 1;
    
    % 0 = dont show figures, but save them to pdf
    % 1 = show figures, and save pdf
    % 2 = show figures, and save eps annotated with parameters
    plotFigures = 1;

    % Do you want to save files to HDF5 (larger file size, but open
    % standard) or use the default matlab format.
    HDF5 = false;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Parameters used for the Koulakov model

    % Activity parameters
    aAct = 0.03;            % Compare, N=100 in Tripplet, their a = 0.03*100
                            % (also used by updateNeighboursTables)
    bAct = 0.11;
    
    gammaAct = 0.05/2; 
    actScaling = NaN;  % This should be scaled by a factor 1e4/nSC * 1e4/nRGC
    
    % Competition parameters
    AComp = -500; % -500 in Tripplet % Note sign
    BComp = 1;
    DComp = 1;
    EComp = 0; % 0 in Tripplet
    alphaComp = 0.5;
    betaComp = 2;
    deltaComp = 2;
    epsilonComp = 0.5;
    
    % Should we mask gradients?
    kMask = 0; % 0 = no, 1 = use reaction scheme,
               % 2 = use subtractive masking
    
    % Chemical parameters
    % !!! Different formulation from Tsigankov 2006 paper

    typeFlag = 1; % 1 = Forward only, 2 = Forward/Backward, 3 = Servo
                  % 4 = reaction
    
    % If you are using typeFlag = 4 you need to recalculate the
    % lookup table whenever you change the alpha or beta values
    
    alphaForwardChem = 20;
    betaForwardChem = 30;

    alphaReverseChem = 20; % !!! Reverse signaling should be weaker !!!
    betaReverseChem = 30;

    alphaServoChem = 20;
    betaServoChem = 30;

    % How much poisson noise in Eph and ephrin systems
    RGCnoiseLevelN = 20;
    SCnoiseLevelN = 20;
    
    % Lookup table for exact chemical reactions
    
    lookupFileChem = 'reactionLookup.mat';

    gridMaxChem = 2;   % Eph and ephrin are in range 0-1, but noise
                       % can make the values larger
    
    nGridChem = 50;    % Number of points per dim in the 4D lookup
                       % table, 20 gives 1e-3 abs error, 
                       % 50 gives 1e-4 abs error
    
    
    KtfChem = 10;    % Trans forward
    KtrChem = 10;    % Trans reverse
    KcaChem = 5;     % Cis axon
    KcdChem = 5;     % Cis dendrite
    
    chemInteractionLookup = []; % RGC are columns, SC are rows
    % It contains the total chemical component
    % alpha*forwardReaction + beta*forwardReaction
    % + alpha*reverseReaction + beta*reverseReaction
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Retinal information
    
    eyeType = 'disk'; % sphere or disk
    RGCdensity = 'uniform'; % 'uniform' or 'varying'
    
    gradientGenerationMethod = 'naive'; % phenotype, diffusion or naive
    phenotype = 'WT';
    
    % Location of RGC in the retina and in SC
    RGCnt = []; % Nasal-Temporal retinal coordinate: 0-1 (A)
    RGCdv = []; % Dorsal-Ventral retinal coordinate: 0-1 (B)

    RGCntPadding = [];
    RGCdvPadding = [];
    
    % Spherical Retinal RGC coordinates
    RGCphi = [];
    RGCtheta = []; % Theta = pi at center of eye opening

    paddingPhi = []; % Keep the rejected positions, useful for triangulation
    paddingTheta = [];

    maxTheta = pi/180*(90+37); % Dan data P0 (Sterratt mail 2 dec 2011)

    RGCEphA = []; % Forward signaling
    RGCEphB = [];
    RGCephrinA = []; % Reverse signaling
    RGCephrinB = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Parameters for modelling different RGC subtypes, with higher
    % correlation within than between the population
    
    nRGCtypes = 1;
    RGCtype = []; % cell type, subpopulation ID
    RGCdPop = 0.1; % Distance at which within subpopulation
                    % correlation falls below between subpopulation
                    % peak correlation.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Parameters relevant when we only allow local jumps,
    % i.e. synapses can only be created near existing synapses
    % belonging to the same axons.
    
    % Where the RGC axon enters SC
    RGCml = []; % Medial-Lateral SC coordinates: 0-1 (B)
    RGCwidth = 2; % How wide is the RGC arbor (in ml direction)
                  % OBS: fraction of total width (scaled)
                  % Value must be 2 or larger to guarantee
                  % behaviour similar to Koulakov if stepFast.c is
                  % used in stepMex.m

    RGCmlSpread = 0.25;
    
    useLocalJumps = false; % updateNeighboursTable is affected
    
    maxSynapseJumpLength = 0.1; % 0.1 lets entire map form in 10000
                                % iterations without acitivyt cues
                                % at 0.2 the local jumps are not
                                % affecting final termination zone
                                % size for case without activity
                                % and chemical cues.
    
    
    nInitSynapses = 20; % Only used if useLocalJumps = true
                        % The synapses will be placed within
                        % a segment of RGCwidth around axon
    
    % If local jumps are used, ie stepFastLocal.c, please set
    % RGCwidth = 0 to conserve memory. It hogs a lot of resources!
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % SC information
    
    % Location of neurons in SC
    SCap = []; % Anterior-Posterior SC coordinates: 0-1 (A)
    SCml = []; % Medial-Lateral SC coordinates: 0-1 (B)
    SCapPadding = [];
    SCmlPadding = [];
    
    SCz = []; % Z coordinate in SC
    SCzPadding = [];
    
    SCdepth = 0.4;
    SC3Dflag = false;
    
    SCscale = 2.072; % Scale to turn normalized units into mm
                     % WT 60 days, Dan Lyngholm, personal
                     % communication, 18 july 2012

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    SCephrinA = []; % Forward signaling
    SCephrinB = [];
    SCEphA = []; % Reverse signaling
    SCEphB = [];

    Isl2PositiveRGC = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Neighbour relations in the SC
    neighbourSC = {};  % For each SC, lists nearby SC neurons
                       % (including self)
    nNeighbourSC = [];
    
    SCmask = []; % Used by filterSCPosition.m
    
    synapseNeighbourhood = {}; % Used for local jumps, neighbours
                               % around SC neuron that a new RGC
                               % axon synapsing on this particular
                               % neuron can jump to.
    
    nSynapseNeighbourhood = [];

    % Synapse connectivity
    presynapticConnections = [];  % One column per SC, NaN padded
    
    presynapticWeight = [];       % One column per SC, weight of each connection
                                  % W=3 means a pair of neurons has
                                  % three synapses this speeds things up.
    
    %%%% The stepFast.c framework does not use these three variables below
    
    postsynapticConnections = []; % One column per RGC, NaN padded
    postsynapticWeight = [];      % One column per RGC, weight of
                                  % each connection
    
    connectionMatrix = [];        % One column per RGC, one row per SC
    
    %%%%
         
    % How many connected RGC axons (per SC neuron)
    numPresynapticConnections = []; 

    % How many connected SC neurons (per RGC neuron)
    numPostsynapticConnections = [];
    
    % How many target SC synapses for each RGC axon
    totalWeightRGC = []; 

    % Total incomming weight for each SC
    totalWeightSC = [];         

    % How many different RGC can each SC be connected to at most
    maxConnections = 100;    % NaN = auto, allows connections to
                             % all neighbours within range
                             % auto only works for global jumps

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Activity look up tables
    CAct = [];
    UAct = {};

    % Servo precalc
    servoExp = exp(-1); % !!! This one has to be modified depending
                        % on gradients.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % At the actionIteration the actionHandler is called, both are vectors
    
    actionHandler = {};
    actionIteration = [];
    actionParameters = {};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    RGCalive = [];
    
    RGCcolour = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    info = []; % Empty variable to store information about the simulation
    
  end

  methods
 
    % Most functions are defined in separate file



   	% Constructor and destructor defined in classdef file (ie here)

    function obj = RetinalMap(configFile)

      if(exist('configFile') & length(configFile) > 4 & strcmpi(configFile(end-3:end),'.mat'))
         % The config file was a mat file, load the state instead
         obj.loadState(configFile);
         return
      end
      
      disp('Constructor called')
      
      obj.initializeRandomGenerator();

      obj.simID = ceil(rand(1)*1e12);
      obj.timeStamp = now();
      
      if(exist('configFile'))
        obj.loadExperimentConfigFile(configFile);
      else
        disp('No config file specified, no initialisation done.')
        return
      end
      
      % Make sure the directories needed exist, otherwise create them
      obj.checkDirectories();
      
      fprintf('Using %d RGC and %d SC\n', obj.nRGC, obj.nSC)

      obj.placeNeurons();
      obj.setupGradients();

      if(isnan(obj.actScaling))
        obj.actScaling = 1e4/obj.nSC;
        fprintf('Scaling external correlations by a factor %d\n', obj.actScaling)
      else
        fprintf('User overriding external correlation scaling, factor %d\n', obj.actScaling)
      end
           
      if(strcmp(obj.eyeType,'disk') & obj.plotFigures)
        obj.plotGradients1D();      
      end
        
      if(obj.typeFlag == 4)
        % We need to precalculate the tables
        obj.makeReactionLookupTable();
      end
        
      % Set up initial connections, and neighbourhood lookup tables
      obj.setupInitialConnectivity();

     
      if(obj.plotFigures)
        obj.plotMapForward();
      end
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function delete(obj)
      disp('Destructor called')
      % close all
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % All other functions, defined in separate files

    %%% Initialisation functions

    initializeRandomGenerator(obj)
    loadExperimentConfigFile(obj,configFile)    

    placeRetinalRGCdisk(obj,minNum)
    [NT,DV,NTpad,DVpad,mask] = filterRetinalPositionsDisk(obj,NT,DV)
    setRGCdensity(obj)

    placeRGCdiskRealDensity(obj,n)    
    placeRetinalRGCdensity(obj,nRGC)

    setNumberOfRGC(obj,nRGC,idx)
    setNumberOfSC(obj,nSC)

    hyperSphere(obj,nPoints)
    placeRetinalRGCsphere(obj,minNum)
    [theta,phi,thetaReject,phiReject,mask] = ...
      filterRetinalPositionSphere(obj,thetaRaw,phiRaw)

    placeSC(obj,make3D)
    placeSCregularGrid(obj)
    [AP,ML,APpad,MLpad,mask] = filterSCPosition(obj,AP,ML,cropFlag);
    [AP,ML,Z,APpad,MLpad,Zpad,mask] = filterSCPosition3D(obj,AP,ML,Z)

    makeNaiveGradients(obj, varargin)
    maskGradients(obj, plotFlag)
    maskGradientsSubtractive(obj)
    
    setupGradientsRetina(obj)
    setupGradientsSC(obj)

    loadGradients(obj,phenoType); % Loads gradients from DS

    steadyStateConc = calculateGradient(obj, ...
                                        edges, ...
                                        edgeLen, ...
                                        crossSection, ...
                                        nodeVolume, ...
                                        pointSources, ...
                                        pointSourceValues, ...
                                        diffusionConstant)

    signalNoise = addPoissonNoise(obj,signal,N)
    
    setupInitialConnectivity(obj)
    randomizeRGCAxonPositionWithBias(obj)
    randomizeRGCAxonPositionNoBias(obj)

    placeInitialSynapses(obj)

    % Allows you to cut away part of the retina or SC, or rotate
    % part of the retina - not yet adding two nasal parts.
    surgery(obj,type)
    
    
    %%% Helper functions

    checkDirectories(obj)
    
    [X,Y,Xout,Yout] = dMin(obj,constrainFcn,nPoints,dMin,packingFactor)    
    [X,Y,Z,Xout,Yout,Zout] = dMin3D(obj,constrainFcn,nPoints, dMin,packingFactor)
    
    d = retinalDistance(obj) % These are on the sphere
    d = SCDistance(obj)
    
    [edges, edgeLen, crossSection, nodeArea] = ...
       getTriangulationRetinalSphere(obj)
    [edges, edgeLen, crossSection, nodeArea] = ...
       getTriangulationSC(obj)
    [edges, edgeLen, crossSection, nodeArea] = ...
       getTriangulation(obj,x,y,z,xPad,yPad,zPad)
    
    % This version of the function is used if we want to clamp the
    % edges to a value, because then the padding needs to be
    % included in the calculation. Normally we exclude it.
    [edges, edgeLen, crossSection, nodeArea] = ...
    getTriangulationRetinalSphereIncludingPadding(obj);

    [edges,edgeLen] = getAllEdges(obj,x,y,z)
    crossSection = getAllCrossSections(obj,x,y,z,convhull,edges)

    % X = theta, Y = phi; or X = NT, Y = DV, only use subsetIdx if
    % you want to restrict which neurons you can possibly select
    closestIdx = findClosestRetinalNeuron(obj, X, Y, subsetIdx)
    neuronIdx = findRetinalNeurons(obj,theta,phi,radius)
    closestIdx = findSCneuron(obj,ml,ap,d)

    [preIdx,weight] = findPresynapticRGC(obj, SCidx)    
    
    [x,y,z] = getCartesianCoordinatesRGC(obj)
    
    updateNeighboursTables(obj)
    calculateCU(obj)
    calculateCU3D(obj) % Used for 3D SC

    modifyActivityComponent(obj,gammaAct,bAct,aAct)

    % Calculate a sharper U drop, with a bigger basin
    calculateAlternativeU(obj)

    
    % Used to calculate the chemical reaction exactly
    makeReactionLookupFile(obj, matlabPoolSize)
    makeReactionLookupTable(obj)

    [RaLd, LaRd,RaLa, RdLd] = ...
        receptorLigandReactionSteadyState(obj, Ra0, La0, Rd0, Ld0)


    RGCidx = getRandomNewRGC(obj,SCidx)
    RGCidx = getRandomExistingConnection(obj,SCidx)
    addSynapse(obj, SCidx, RGCidx)
    removeSynapse(obj, SCidx, RGCidx)
    
    [Etotal,Echem,Eact,Ecomp] = calculateTotalEnergy(obj)
    [Eadd,Echem,Eact,Ecomp] = calculateSynapseAdditionEnergy(obj,SCidx,RGCidx)
    [Eremove,Echem,Eact,Ecomp] = calculateSynapseRemovalEnergy(obj,SCidx,RGCidx)
    Esynapse = calculateSynapseActivityEnergy(obj,SCidx,RGCidx,addFlag)
    Esynapse = calculateSynapseCompetitionEnergy(obj,SCidx,RGCidx,addFlag)
    Esynapse = calculateSynapseChemicalEnergy(obj,SCidx,RGCidx)

    saveState(obj, fileName)
    loadState(obj, fileName, skipTables)
    saveIter(obj, fileName)
    iterSequence = loadIterSequence(obj,skipTables,iter,dataPath,fileNameMask,uniqueFlag)
    seqFilt = filterSequence(obj,seq,stepSize,skipZeroFlag,keepFirstNonZeroFlag)
    
    [RGCcentAP,RGCcentML] = RGCprojectionCentroids(obj)
    [RGCmaxAP,RGCmaxML] = RGCprojectionMax(obj,useTableFlag)
    [SCmaxNT,SCmaxDV] = SCprojectionMax(obj,useTableFlag)

    col = getSphereColour(obj,theta,phi)

    [RGCcol,SCcol] = getNeuronColoursFromImage(obj,refImageName,flipImg)
    
    % Test functions using the global energy instead
    [Etotal,Echem,Eact,Ecomp] = calculateTotalEnergyGlobal(obj)

    % Export and import functions
    importConnectionMatrix(obj,connectionFile,stateFile)
    importConnectionIndexes(obj,indexFile,weightFile,stateFile)
    exportGradients(obj)
    exportMapToDavid(obj);
    
    % This function converts between a presynaptic or postsynaptic
    % representation of the connectivity
    convertConnectionTables(obj,action);
    
    RGCcolour = getRGCcolour(obj)
    
    [k, kLow, kHigh] = fitErrorEstimate(obj,x,y,func,k0,nRep)
    
    % yLow is 2.5% percentile, and yHigh is 97.5% percentile
    [xRange,yMedian,yLow,yHigh,yAll] = fitErrorEstimateN(obj,x,y,func,k0,nRep,nPoints)   

    [k,kAll,xRange,yMedian,yAll] = fitJackKnife(obj,x,y,func,k0)    
    [xRange,yMedian,yAll] = loessJackKnife(obj,x,y)

    shareAxis(figs)
        
    %%% Run functions
    
    step(obj,nSteps)
    [stepCtr, stepEnergy] = testconvergence(obj)

    % Wrapper function to Mex function
    stepMex(obj,nSteps)
    
    % Mex function, call it through the stepMex wrapper function
    % (global jumps)
    [presynapticConnections, presynapticWeight, ...
     numPresynapticConnections, totalWeightRGC, ...
     totalWeightSC, time] = stepFast(obj,nSteps);    

    % Mex function, call it through the stepMex wrapper function
    % (local jumps)   
    [presynapticConnections, presynapticWeight, ...
     numPresynapticConnections, totalWeightRGC, ...
     totalWeightSC, time] = stepFastLocal(obj,nSteps);    

    % Implements the reaction scheme + local jumps
    [presynapticConnections, presynapticWeight, ...
     numPresynapticConnections, totalWeightRGC, ...
     totalWeightSC, time] = stepFastLocalChem(obj,nSteps);    
    
    % Experimental Data
    [WT,B2KO,SCsize] = getDanSegregationFits(obj)
    
    % Ation handler that takes input
    % aParam = {{'varName1',varValue1},{'varName2',varValue2},...}
    % and sets obj.varName1 = varValue1; etc.
    actionHandlerSetParam(obj,aParam);
    
    %%% Analysis functions
    
    allEqual = paramDiff(obj,otherObj)
    
    plotRetinaRGCdisk(obj)
    plotRetinalRGCsphere(obj,refNeuron)
    fig = plotGradients1D(obj)
    fig = plotGradients2D(obj, invertSCephrinA, forwardOnly, markIsl2)
    fig = plotEphrinInSC(obj)
    
    plotMesh(obj, x, y, varargin)

    plotRetinalGradient3D(obj, gradient, color)
    plotRetinalGradient2D(obj)
    plotRetinalDistVsGradient(obj, origoNeuron, gradient)
    plotRetinalGradientDisk(obj)
    
    plotMap(obj)
    [RGCidx,SCidx] = makeVirtualInjection(obj,SCap,SCml,injRadius) 
    [fig,sp] = plotVirtualInjection(obj,SCidx,RGCidx,SCidx2,RGCidx2)

    
    fig = plotMapReverse(obj,interactiveFlag,SCcol,RGCcol,markIsl2)
    fig = plotMapForward(obj,interactiveFlag,markIsl2)
    fig = plotMapImage(obj,refImageName,interactiveMap,flipImg)

    
    distHist = mapPairDistanceHistogram(obj,nBins)    
    fig = plotAxisProjection(obj, showAxis, useAll, reuseFigFlag,plotType,saveFigFlag,weightScale,alpha,showExpFlag)
    plotAxisMap(obj, axisDir, nBins, axisRange, colorbarFlag, plotFlag);
    plotBothAxisMaps(obj)
    plotAxonTerminationCentroids(obj,onlyPlotSCflag)
    fig = plotTerminationZoneSize(obj,sequence,percentile,titleText)
    fig = plotTerminationZoneSizeCompare(obj,sequences,percentile,legendText,colours)
    plotSC(obj,noAxisFlag)
    plotSCinputIdentity(obj)
    plotInputClustering(obj)

    % Pass kEst along if it is already known
    [areaCoverage,kEst,fig] = ...
      contourAnalysis(obj,RGCidx,percentile,RGCweights,SCidx,plotFlag,kEst,approxFlag)

    [areaCoverage,kEst,fig] = contourAnalysisKDE(obj, RGCidx, ...
                                                 percentile, plotFig)
    
    makeMapMovie(obj,iterSequence, RGCidx, SCidx, resolution, frameRate)
    
    fig = plotNumberOfSynapses(obj,plotType,oldFig)
    inspectPotential(obj,addFlag,SCidx)
    inspectSelfAddCost(obj)

    [nearestNeighbourRegularityIndex, voronoiRegularityIndex] = ...
        nearestNeighbourStatistics(obj,structure,X,Y)
 
    [distance, segregation,k,kLow,kHigh] = ...
        virtualInjectionSegregationExperiment(obj, injectionType, nRep, maxDist, k0);
    [dEadd,dEremove] = sampleRelativeEnergyStrength(obj, nSamples, type);
    [hasEctopic,fig] = virtualInjectionEctopicExperiment(obj,RGCnt, RGCdv, injectionRadius);    
    fig = virtualInjectionMappingExperiment(obj,axis)
    [SCdistance,retinaFlipFlag] = virtualInjectionVectorExperiment(obj,axisDir,nRep) 
    fig = plotVirtualExperimentVector(obj,axisDir, titleStr)

    NT = locateCollapsePoint(obj,debugFlag)    
    
    [MCscore, MCratio] = calculateMapCoherency(obj)
    fracCovered = calculateSCcoverage(obj,synapseFrac)
    
    [TZspreadInj,nRGCLabeled,nSCLabeled] = ...
        calculateTerminationZoneSizeInjection(obj,nRep,percentile,plotFlag)

    plotMapSpread(obj, interactiveFlag)
    
    [grid,fig] = makeProjectionGrid(obj, nPoints, radius, debugFigs)
    [H,majorAng,mmRatio,fig] = superposedProjection(obj,RGCidx,showRange)
    [dFuzz,kFuzz,cFuzz95,fig] = plotMapFuzzyness(obj, axisDir, nRep);
    dFuzz = calculateMapFuzzynessDistance(obj,axisDir)
    [k,c] = calculateMapFuzzynessRegion(obj,axisDir)

    [k,fig] = plotVirtualExperimentSegregation(obj,axisDir,titleText, ...
                                                    refCurves,refColour,lengthScale,legendText)
    
    [k,iter,fig] = ...
        plotVirtualExperimentSegregationSequence(obj,axisDir,sequence,plotPoints,titleText)

    placeIterationsOnDevelopmentalTimeline(obj)
    
    matchIterationToAgeSegregationExperiment(obj,injectionType,nRep)
    [age,bestIter] = matchIterationToAgeSegregationExperimentScore(obj)
    [age,bestIter] =  matchIterationToAgeSegregationExperimentDist(obj,normaliseSize)    
    
    % Used for supplementary figure, B2KO article
    [bestIter,bestIdx,kFit,kMatch] = findSegregationCrossing(obj,matchObj,kMatch,stopAtIter)
    
    % Analyse seggregation data from Dan
    [age,kWT,kB2KO, kWTJK, kB2KOJK] = fitSegregationData(obj,scaleDist,visibleFlag)
    [k,kLow,kHigh] = fitSegregation(obj,distance,segregation,k0)
    
    checkTopology(obj, rMin, rMax, RGCidx)

    % Display the parameters under a figure, assuming normal layout
    addParametersToFigure(obj,fileName,textString)
    fig = plotParameters(obj,fontsize)
    fig = plotQuantificationPanel(obj, report,fontsize)
    
    trimFigure(obj,fig)
    
    % Copy a set of figures over to subplots in one figure
    fig = groupPlots(obj,subplotY,subplotX,figHandleList,singleFlag, ...
                     fontsize,markersize,markerSizeDot,linewidth, ...
                     subplotPositionIdx)
    
    saveMultiPagePS(obj,fileName,fig,appendFlag,format)
    
    report = makeReport(obj, plotFigures, DWgrid)
    
    loadAndDo(obj,directory,pattern,action);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Some new test functions
    
    stepNew(obj,nSteps,nSynPerStep,nPatternsPerStep,waveSize,dt)
    
  end

end
