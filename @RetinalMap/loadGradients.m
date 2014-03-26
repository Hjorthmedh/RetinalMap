function loadGradients(obj,phenoType)

  % This loads the table that David Sterrat had prepared
  % Added the B-system
  % infoFile = 'gradients/Eph-ephrins-both-axis-JH.csv';

  fid = fopen(obj.gradientInfoFile,'r');
  fprintf('Reading gradient info file: %s\n', obj.gradientInfoFile)
  
  header = fgets(fid);
  
  if(header(end) == 13)
    % This works around a bug in textscan in matlab 2013a
    header = header(1:end-1);
  end
  
  columnTokens = textscan(header,'%s','delimiter',',');
  columnTokens = strrep(columnTokens{1},'"','');
  
  % Clear the old gradients
  obj.RGCEphA = zeros(obj.nRGC,1);
  obj.RGCEphB = zeros(obj.nRGC,1);
  obj.RGCephrinA = zeros(obj.nRGC,1);
  obj.RGCephrinB = zeros(obj.nRGC,1);
  
  obj.SCephrinA = zeros(obj.nSC,1);
  obj.SCephrinB = zeros(obj.nSC,1);
  obj.SCEphA = zeros(obj.nSC,1);  
  obj.SCEphB = zeros(obj.nSC,1);    
  
  RGCEphAMolecule = {};
  RGCEphBMolecule = {};
  RGCephrinAMolecule = {};
  RGCephrinBMolecule = {};

  SCEphAMolecule = {}; 
  SCEphBMolecule = {};
  SCephrinAMolecule = {};
  SCephrinBMolecule = {};

  RGCEphAPartial = [];
  RGCEphBPartial = [];
  RGCephrinAPartial = [];
  RGCephrinBPartial = [];

  SCEphAPartial = [];
  SCEphBPartial = [];
  SCephrinAPartial = [];
  SCephrinBPartial = [];
  
  RGCEphAcol = [];
  RGCEphBcol = [];
  RGCephrinAcol = [];
  RGCephrinBcol = [];
  
  SCEphAcol = [];
  SCEphBcol = [];
  SCephrinAcol = [];
  SCephrinBcol = [];
  
  
  % Now add the components as we go through the file
  
  infoStr = fgets(fid);
  while(infoStr ~= -1)
  
    % First, check if this gradient should be included
    weight = getParameter(infoStr,phenoType);
    present = getParameter(infoStr,'Present');

    if(weight & present)
      % Ok, we need this gradient
      
      C0 = getParameter(infoStr,'C0');
      C1 = getParameter(infoStr,'C1');
      C2 = getParameter(infoStr,'C2');
      C3 = getParameter(infoStr,'C3');      
      
      axisDir = getParameter(infoStr,'Axis');
      
      switch(axisDir)
        case 'NT'
          X = obj.RGCnt;
        case 'DV'          
          X = obj.RGCdv;          
        case 'AP'
          X = obj.SCap;
        case 'ML'
          X = obj.SCml;
        otherwise
          fprintf('Unknown axis: %s\n', axisDir)
          keyboard
      end

      try
        val = C0 + C1*exp(-C2*abs(X-C3));

        % Make sure it is non-negative -- this gives us more
        % freedome when parameterising
        
        val = max(0,val);
        
        if(weight < 1)
          mask = rand(size(val)) < weight;
          val = val .* (mask);
        end
        
      catch e
        getReport(e)
        keyboard
      end
        
      % Now we need to add it to the relevant gradient
      structure = getParameter(infoStr,'Structure');
      molecule = getParameter(infoStr,'molecule');
      
      switch(structure)
        
        case 'Retina'

          if(strcmpi(molecule(1:4),'EphA'))
            disp('Adding EphA to retina')
            obj.RGCEphA = obj.RGCEphA + val;
            
            RGCEphAMolecule{end+1} = molecule;
            RGCEphAPartial(end+1,:) = val;
            RGCEphAcol(end+1,:) = getParameter(infoStr,'plotColour');

            % Check if it was Isl2 that were added
            if(strcmpi('EphA3-Isl2',molecule))
              disp('Saving indexes of Isl2 positive RGC cells')
              obj.Isl2PositiveRGC = find(val);
            end
            
          elseif(strcmpi(molecule(1:4),'EphB'))
            disp('Adding EphB to retina')
            obj.RGCEphB = obj.RGCEphB + val;
            
            RGCEphBMolecule{end+1} = molecule;
            RGCEphBPartial(end+1,:) = val;
            RGCEphBcol(end+1,:) = getParameter(infoStr,'plotColour');
            
          elseif(strcmpi(molecule(1:7),'ephrinA'))
            disp('Adding ephrinA to retina')
            obj.RGCephrinA = obj.RGCephrinA + val;

            RGCephrinAMolecule{end+1} = molecule;
            RGCephrinAPartial(end+1,:) = val;
            RGCephrinAcol(end+1,:) = getParameter(infoStr,'plotColour');
            
          elseif(strcmpi(molecule(1:7),'ephrinB'))
            disp('Adding ephrinB to retina')
            obj.RGCephrinB = obj.RGCephrinB + val;
            
            RGCephrinBMolecule{end+1} = molecule;
            RGCephrinBPartial(end+1,:) = val;
            RGCephrinBcol(end+1,:) = getParameter(infoStr,'plotColour');
            
          else
            fprintf('Unknown molecule %s.\n', molecule)
            assert(0)
          end
          
        case 'SC'

          if(strcmpi(molecule(1:4),'EphA'))
            disp('Adding EphA to SC')
            obj.SCEphA = obj.SCEphA + val;

            SCEphAMolecule{end+1} = molecule;
            SCEphAPartial(end+1,:) = val;
            SCEphAcol(end+1,:) = getParameter(infoStr,'plotColour');
            
          elseif(strcmpi(molecule(1:4),'EphB'))
            disp('Adding EphB to SC')
            obj.SCEphB = obj.SCEphB + val;

            SCEphBMolecule{end+1} = molecule;
            SCEphBPartial(end+1,:) = val; 
            SCEphBcol(end+1,:) = getParameter(infoStr,'plotColour');           
            
          elseif(strcmpi(molecule(1:7),'ephrinA'))
            disp('Adding ephrinA to SC')
            obj.SCephrinA = obj.SCephrinA + val;
            
            SCephrinAMolecule{end+1} = molecule;
            SCephrinAPartial(end+1,:) = val;
            SCephrinAcol(end+1,:) = getParameter(infoStr, 'plotColour');
            
          elseif(strcmpi(molecule(1:7),'ephrinB'))
            disp('Adding ephrinB to SC')
            obj.SCephrinB = obj.SCephrinB + val;

            SCephrinBMolecule{end+1} = molecule;
            SCephrinBPartial(end+1,:) = val;            
            SCephrinBcol(end+1,:) = getParameter(infoStr, 'plotColour');
          else
            fprintf('Unknown molecule %s.\n', molecule)
            assert(false)
          end
          
        otherwise
          
          fprintf('Unknown structure: %s\n', structure)
          assert(false);
          
      end
      
    end
    
    infoStr = fgets(fid);    
  end
  

  % Normalise the signal so peak is at 1 for WT (Reber)
  obj.RGCEphA = obj.RGCEphA / 3.54;
  % obj.RGCephrinA = obj.RGCephrinA / 1;
  % obj.RGCEphB = obj.RGCEphB / 1;
  % obj.RGCephrinB = obj.RGCephrinB / 1; 

  % obj.SCephrinA = obj.SCephrinA / 1;
  % obj.SCEphA = obj.SCEphA / 1;
  % obj.SCephrinB = obj.SCephrinB / 1;
  % obj.SCEphB = obj.SCEphB / 1;
  

  % Normalise the partial plots also
    
  RGCEphAPartial = RGCEphAPartial / 3.54;
  % RGCephrinAPartial = RGCephrinAPartial / 1;
  % RGCEphBPartial = RGCEphBPartial / 1;
  % RGCephrinBPartial = RGCephrinBPartial / 1;
    
  % SCephrinAPartial = SCephrinAPartial / 1;  
  % SCEphAPartial = SCEphAPartial / 1;
  % SCephrinBPartial = SCephrinBPartial / 1;
  % SCEphBPartial = SCEphBPartial / 1;
  
  if(obj.RGCnoiseLevelN | obj.SCnoiseLevelN)
    fprintf('Adding noise, noise level = %f (RGC) %f (SC)\n', ...
            obj.RGCnoiseLevelN, obj.SCnoiseLevelN)
  end
  
  if(obj.RGCnoiseLevelN)
    obj.RGCEphA = obj.addPoissonNoise(obj.RGCEphA,obj.RGCnoiseLevelN);
    obj.RGCEphB = obj.addPoissonNoise(obj.RGCEphB,obj.RGCnoiseLevelN);
    obj.RGCephrinA = obj.addPoissonNoise(obj.RGCephrinA, ...
                                         obj.RGCnoiseLevelN);
    obj.RGCephrinB = obj.addPoissonNoise(obj.RGCephrinB, ...
                                         obj.RGCnoiseLevelN);
  end
          
  if(obj.SCnoiseLevelN)
    obj.SCephrinA = obj.addPoissonNoise(obj.SCephrinA, ...
                                        obj.SCnoiseLevelN);
    obj.SCephrinB = obj.addPoissonNoise(obj.SCephrinB, ...
                                        obj.SCnoiseLevelN);  
    obj.SCEphA = obj.addPoissonNoise(obj.SCEphA,obj.SCnoiseLevelN);  
    obj.SCEphB = obj.addPoissonNoise(obj.SCEphB, ...
                                     obj.SCnoiseLevelN);
  end

  
  if(obj.plotFigures)
    
    if(~exist(obj.figurePath))
      mkdir(obj.figurePath);
    end
    
    coord = { obj.RGCnt, obj.RGCdv, ...
              obj.RGCnt, obj.RGCdv };
    
    xLabelName = { 'Nasal-Temporal', 'Dorsal-Ventral', ...
                   'Nasal-Temporal', 'Dorsal-Ventral' };
    gradient = { RGCEphAPartial, RGCEphBPartial, ...
                 RGCephrinAPartial, RGCephrinBPartial};
    legText = { RGCEphAMolecule, RGCEphBMolecule, ...
                RGCephrinAMolecule, RGCephrinBMolecule };

    plotColour = { RGCEphAcol, RGCEphBcol, RGCephrinAcol, RGCephrinBcol };
    
    makeAreaPlot(coord, gradient, xLabelName, legText, plotColour);
    fName = sprintf('%s/%s-retinalGradients.pdf', obj.figurePath, obj.simName);
    saveas(gcf,fName,'pdf')
    
    xLabelName = { 'Anterior-Posterior', 'Medial-Lateral', ...
                   'Anterior-Posterior', 'Medial-Lateral' };
    
    coord = { obj.SCap, obj.SCml, ...
              obj.SCap, obj.SCml };
    
    gradient = { SCephrinAPartial, SCephrinBPartial, ...
                 SCEphAPartial, SCEphBPartial };
    legText = { SCephrinAMolecule, SCephrinBMolecule, ...
                SCEphAMolecule, SCEphBMolecule };

    plotColour = { SCephrinAcol, SCephrinBcol, SCEphAcol, SCEphBcol };    
    
    makeAreaPlot(coord, gradient, xLabelName, legText, plotColour);    
    fName = sprintf('%s/%s-SCGradients.pdf', obj.figurePath, obj.simName);

    switch(obj.plotFigures)
      case {0,1}
        saveas(gcf,fName,'pdf');
      case 2
        obj.addParametersToFigure(strrep(fName,'.pdf','.eps'));
    end  
    
  end
    
  obj.plotGradients1D();

  fclose(fid);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function d = getParameter(str,columnName);

    % Find the column
    c = find(strcmpi(columnTokens,columnName));

    if(isempty(c))
      fprintf('Unable to find column: %s\n', columnName)
      disp('Available columns: ')
      for i = 1:numel(columnTokens)
        fprintf('%s ', columnTokens{i})
      end
      fprintf('\n')
      disp('Debug mode, type dbquit to exit')
      keyboard
    end
    
    % Extract the relevant token
    str = strrep(str,',,',', ,'); % strtok doesnt handle emptyness well
    
    for i=1:c;
      [d,str] = strtok(str,',');
    end
    
    d = strrep(d,'"','');
    
    dn = str2num(d);
    

    if(~isempty(dn))
      d = dn;
    end
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function makeAreaPlot(coord, gradient, xLabelName, legText, plotColour)
    
    figure

    for i = 1:4
      
      [~,idx] = sort(coord{i});
      
      subplot(2,2,i)

      
      if(~isempty(gradient{i}))
        p = area(repmat(coord{i}(idx),1,numel(legText{i})), ...
                 transpose(gradient{i}(:,idx)));

        % Here we assume coords are 0 to 1
        if(strcmpi(obj.phenotype,'WT'))
          axis([0 1 0 1]);
        end
        
        legend(p(end:-1:1),legText{i}{end:-1:1},'location','best');
      else
        p = [];
      end
      
      for j = 1:numel(p)
        try
          set(p(j),'facecolor',plotColour{i}(j,:))
        catch e
          getReport(e)
          keyboard
        end
      end
      
      xlabel(xLabelName{i},'fontsize',20);
      ylabel('Concentration','fontsize',20);
      set(gca,'fontsize',14)
            
      box off
     
      %title(sprintf('Max: %f', max(sum(gradient{i}))))
      
    end

  end
  
end