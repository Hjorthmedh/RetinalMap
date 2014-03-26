% This file is the standard Marker Induction simulation with 500
% RGCs and SC cells
function exitcode =  markerInductionBatch(configFile, method, plotFigures)

if (~exist('method')) 
  method = 'C';
end
if (~exist('plotFigures')) 
  plotFigures = 0;
end


% Load Config file
rS = RetinalMap();
rS.loadExperimentConfigFile(configFile);
rS.plotFigures = false;
rS.placeRetinalRGCdisk();
rS.placeSC();
rS.loadGradients(rS.phenotype);

% Time step for Marker Induction
rS.info.dt = 0.1;
%rS.nSteps = 120000/rS.info.dt;
rS.nSteps = 4800/rS.info.dt;
rS.reportStep = 200/rS.info.dt;
rS.plotFigures = plotFigures;

% Constants of Marker Induction
rS.info.alpha = 0.05;
rS.info.beta =  0.01; % /3.64;
rS.info.gamma = 0.1;

% To correspond to our equations, we therefore need from David
% Willshaw's original script:
% kappa=k=MULT*MINRET*RADIUS
rS.info.kappa = 0.72*0.14*1/2;

% Scale up RGC EphA gradients
rS.info.RGCEphAScale = 3.5;

% Set in the ComparePhenotypeModelling script...
if(exist('/disk/scratch/sterratt/projects/rettect/koulakov/'))
  % If on DS system, the path above exists, then use it instead of default.
  rS.dataPath = ['/disk/scratch/sterratt/projects/rettect/koulakov/' rS.dataPath];
end
  
% "Initial pattern of innervation is random" (Fig. 7 caption, p. 2712)
% Initialise random weak matrix
rS.connectionMatrix = 1e-4*rand(rS.nSC, rS.nRGC);
rS.connectionMatrix = rS.connectionMatrix./ ...
    (ones(rS.nSC, 1) * sum(rS.connectionMatrix, 1))
assert(~any(any(isnan(rS.connectionMatrix))))

% Run simulations
rS = runMarkerInduction(rS, method);

% Do connection conversion
connectionMatrix0 = rS.connectionMatrix;
rS.connectionMatrix(rS.connectionMatrix < 0.001) = 0;
rS.convertConnectionTables('mat2pre');
rS.connectionMatrix = connectionMatrix0;
rS.saveState();

exitcode = 0;
