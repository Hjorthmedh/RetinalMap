%
% Phenotype can be WT,Isl2heterozygous,Isl2homozygous,ephrinA2mm, 
% ephrinA5mm, ephrinA2mmA5mm, TKO
% 
% Action can be: run, analyse
%
% expNum is used if you want to run multiple runs with same
% parameter from another script to do statistics. This allows you
% to separate the runs.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function report = ComparePhenotypeModelling(phenotype,action, ...
                                            expNum,plotFigures,kMask,model, varargin)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if(~exist('model'))
    model = 'koulakov';
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if(~exist('expNum'))
    expNum = '';
  else
    if(ischar(expNum))
      % Fix for DS compiled code
      expNum = str2num(expNum);
    end
    expNum = sprintf('-rep-%d', expNum);
  end
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if(~exist('plotFigures'))
    % Show figures!
    plotFigures = 1;
  else
    if(ischar(plotFigures))
      plotFigures = str2num(plotFigures);
    end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if(~exist('kMask'))
    disp('No masking of gradients.')
    kMask = 0; % 0 = no masking, 2 = subtractive masking
  elseif(ischar(kMask))
    % DS compile code fix, he parses the numbers as strings...
    kMask = str2num(kMask);
  end
  
  % Added extra asserts to help DS
  assert(kMask == 0 | kMask == 2 | kMask == 3);
  
  assert(kMask ~= 1); % This is reaction scheme masking, make sure
                      % we dont do this by accident...
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  DWgrid = true; % Do DW grid analysis?
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  if(strcmpi('Math5',phenotype))
    % Do not do DW grid analysis for Math5 - the code crashes
    DWgrid = false;
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  nSC = 2000; %5000;
   
  typeFlag = 1; %2;
  forwardChem = 90; %20; 
  reverseChem = 0; %10; 
  betaMod = 3/2;

  gammaAct = 0.05/8; % Make it weaker (default Koulakov is 0.05/2)  
  
  % See FitCI.m, uses Stafford 2009 data
  % Assums retina is 5mm (flattened out)
  bActB2KO = 1/(5e-3*1260);
  gammaActB2KO = gammaAct*9.275/22.64;  
  
  noiseLevel = 0; % Lower value means higher noise, 20 looks
                  % good! 0 = no noise

  nStep = 10000*nSC;
  repStep = 1000*nSC; % inf

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  baseSimName = 'ComparePhenotype';
  simName = sprintf('%s-%s-%s%s',baseSimName,model,phenotype,expNum);
  
  experimentDir = sprintf('experiments/%s', baseSimName);
  dataPath = sprintf('SAVE/%s', baseSimName);
  figPath = sprintf('FIGS/%s', baseSimName);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  switch(action)
      
    case 'run'
      
      if(~exist(experimentDir))
          mkdir(experimentDir);
      end
      
      if(~exist(figPath))
        mkdir(figPath)
      end
      
      if(~exist(dataPath))
        mkdir(dataPath)
      end
      
      filename = sprintf('%s/%s.txt', experimentDir, simName);
      writeConfig(phenotype,filename);

      switch(model)
        
        case 'template'

          % This sets up the intial conditions
          r = RetinalMap();
          r.loadExperimentConfigFile(filename);
          r.placeRetinalRGCdisk();
          r.placeSC();
          r.loadGradients(r.phenotype);
           
          % This is an example of how to integrate the your own model with
          % the framework
          templatemodel(r);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        case 'koulakov'
          
          r = RetinalMap(filename);
          r.run();
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
        case 'Gierer2D'
          
          % Setup neuron position and gradients
          r = RetinalMap();
          r.loadExperimentConfigFile(filename);
          r.placeRetinalRGCdisk();
          r.placeSC();
          r.loadGradients(r.phenotype);
              
          % Export gradients and neuron locations to file
          r.exportGradients();
      
          % Write config file for Gierer2D R code
          RGCfile = sprintf('%s-retinal-gradients.txt',r.simName);
          SCfile = sprintf('%s-SC-gradients.txt',r.simName);
          paramFilename = sprintf('%s-parameters.R',r.simName);
          paramFilenameFull = sprintf('%s/%s',r.dataPath, paramFilename);
          
          % The R code runs in the SAVE/Gierer2D code, so no need for path
          outputFile = sprintf('%s-weights.dat',r.simName);
          outputFileFull = sprintf('%s/%s',r.dataPath,outputFile);
      
          fid = fopen(paramFilenameFull,'w');
          fprintf(fid,'% Automatically generated by ComparePhenotypeModelling.m\n');
          fprintf(fid,'rgc.file  = "%s"\n', RGCfile);
          fprintf(fid,'sc.file  = "%s"\n', SCfile);
          fprintf(fid,'op.file  = "%s"\n', outputFile);
          
          for i = 1:length(varargin)
            fprintf('Writing to Gierer config file: %s\n', varargin{i})
            fprintf(fid,'%s\n',varargin{i});
          end
          
          fclose(fid);
      
          % Call Gierer2D code in R
          wd = pwd;
          sysCmd = sprintf('cd %s; %s/gierer/inst/runGierer2Dv3.R %s', ...
                           r.dataPath, wd, paramFilename);
          fprintf('Doing: %s\n',sysCmd)
          system(sysCmd)
                     
          % Import resulting weight matrix back into matlab
          % r.importConnectionMatrix(outputFileFull);
          r.importConnectionIndexes(outputFileFull);
      
          r.curStep = NaN;
      
          % Save state, will contain all information at endpoint
          r.saveState();    

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
        case 'WhiteCow'

          % Setup neuron position and gradients
          r = RetinalMap();
          r.loadExperimentConfigFile(filename);
          r.placeRetinalRGCdisk();
          r.placeSC();
          r.loadGradients(r.phenotype);
              
          % Export gradients and neuron locations to file
          r.exportGradients();
      
          % Write config file for Whitelaw and Cowan R code
          % !!! Should the file names include the path?!
          
          RGCfile = sprintf('%s-retinal-gradients.txt',r.simName);
          SCfile = sprintf('%s-SC-gradients.txt',r.simName);
          paramFilename = sprintf('%s-parameters.R',r.simName);
          paramFilenameFull = sprintf('%s/%s',r.dataPath, paramFilename);
          
          outputFile = sprintf('%s-weights.dat',r.simName);
          outputFileFull = sprintf('%s/%s',r.dataPath,outputFile);
      
          fid = fopen(paramFilenameFull,'w');
          fprintf(fid,'% Automatically generated by ComparePhenotypeModelling.m\n');
          fprintf(fid,'rgc.file  = "%s/%s"\n', r.dataPath, RGCfile);
          fprintf(fid,'sc.file  = "%s/%s"\n', r.dataPath, SCfile);
          fprintf(fid,'op.file  = "%s"\n', outputFileFull);
          fclose(fid);
      
          % Call Whitelaw and Cowan code in R
          sysCmd = sprintf('WhitelawCowan/runWhitelawCowan.R %s', ...
                           paramFilenameFull);
          fprintf('Doing: %s\n',sysCmd)
          system(sysCmd)
                     
          % Import resulting weight matrix back into matlab
          r.importConnectionMatrix(outputFileFull);
      
          r.curStep = NaN;
      
          % Save state, will contain all information at endpoint
          r.saveState();    
          
          
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
        case 'Markerinduction'

          r = markerInductionBatch(filename);
          
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        case 'GrimbertCang'
          useAltP = true;
          Pfactor = 1;
          ignoreGradients = false;
          
          r = runGrimbertCang([],useAltP,Pfactor,ignoreGradients,filename);
          
          
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
        otherwise
          disp(['Unsupported model, use: koulakov, Gierer2D, WhitelawCowan,  ' ...
              'GrimbertCang or Markerinduction'])
          report = [];
        
          return
          
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      end
          
      report = r;
              
    case 'analyse'
      
      skipTables = 1;
      
      filename = sprintf('%s/%s.mat', ...
                         dataPath, simName);
      r = RetinalMap();
      r.loadState(filename,skipTables);

      report = r.makeReport(plotFigures,DWgrid);
      
    otherwise
      fprintf('Unknown action: %s\n', action)
      
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function writeConfig(phenotype, filename)
    
    fid = fopen(filename,'w');
        
    fprintf(fid,'% Automatically generated by ComparePhenotypeModelling.m\n');
    
    fprintf(fid,'obj.nSC = %d;\n', nSC);
    
    % !!! Math5
    if(strcmpi('Math5',phenotype))
      fprintf('RGC population reduced to 10%\n')
      fprintf(fid,'obj.nRGC = %d;\n', nSC*0.1);      
    else
      fprintf(fid,'obj.nRGC = %d;\n', nSC);
    end
      
    fprintf(fid,'obj.eyeType = ''disk'';\n');
    fprintf(fid,'obj.gradientGenerationMethod = ''phenotype'';\n');
    fprintf(fid,'obj.phenotype = ''%s'';\n', phenotype);
    
    fprintf(fid,'obj.RGCnoiseLevelN = %f;\n', noiseLevel);
    fprintf(fid,'obj.SCnoiseLevelN = %f;\n', noiseLevel);
    fprintf(fid,'obj.kMask = %f;\n', kMask);
    
    fprintf(fid,'obj.typeFlag = %d;\n', typeFlag);
    fprintf(fid,'obj.alphaForwardChem = %d;\n', forwardChem);
    fprintf(fid,'obj.betaForwardChem = %d;\n', forwardChem*betaMod);
    fprintf(fid,'obj.alphaReverseChem = %d;\n', reverseChem);
    fprintf(fid,'obj.betaReverseChem = %d;\n', reverseChem*betaMod);
    
    if(strcmpi('Beta2KO',phenotype))
      % !!! Beta2 KO
      fprintf(fid,'obj.bAct = %d;\n', bActB2KO);
      fprintf(fid,'obj.gammaAct = %d;\n', gammaActB2KO);          
    else
      % WT
      fprintf(fid,'obj.gammaAct = %d;\n', gammaAct);          
      % Use default bAct value.
    end
    
    fprintf(fid,'obj.simName = ''%s'';\n', simName); 
    fprintf(fid,'obj.dataPath = ''%s'';\n', dataPath);
    fprintf(fid,'obj.figurePath = ''%s'';\n', figPath);    

    fprintf(fid,'obj.nSteps = %d;\n', nStep);

    fprintf(fid,'obj.useLocalJumps = false;\n');      
    fprintf(fid,'obj.RGCwidth = 2;\n');      
      
    % We want intermediate saves    
    fprintf(fid,'obj.reportStep = %d;\n',repStep);
    
    % Display figures
    fprintf(fid,'obj.plotFigures = %d;\n',plotFigures)
    
    % disp('Changing gradient generation file from default.')
    % fprintf(fid,'obj.gradientInfoFile = ''gradients/Eph-ephrins-more-Isl2-paper.csv'';\n')
    
    fclose(fid);
  
  end

  
end
