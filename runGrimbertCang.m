% This code sets up and runs a Grimbert and Cang simulation

function r = runGrimbertCang(phenotype,useAltP,Pfactor,ignoreGradients,experiment)

  
  % Cached data, used by estimateDensity
  Xgrid = [];
  Ygrid = [];
  mask = [];
  densityNormalisation = [];
  SCmaskPad = [];
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if(exist('experiment') & ~isempty(experiment))
    
    % If we use an experiment file, then that specifies the phenotype
    assert(isempty(phenotype));
    
    r = RetinalMap(experiment);

    if(~isempty(phenotype))
      fprintf('Ignoring phenotype option %s, using %s from experiment file.\n', ...
              phenotype, r.phenotype)
    end

    phenotype = r.phenotype;
    
    allowRenameSimName = false;
    
  else
  
    if(~exist('phenotype') | isempty(phenotype))
      phenotype = 'wt';
    end
  
    experiment = sprintf('experiments/GrimbertCang-%s.txt',phenotype);

    allowRenameSimName = true;
    
    r = RetinalMap(experiment);

  end
    

  dt = 1;

  % Arbor placeing parameters
  sigmaE = 0.075;
  sigmae = 0.025;
  sigmaNs = 0.015;
  
  nArborsInit = 7;
  nArbors = 3;
  
  % F-function parameters (Ffunc)
  theta = 0.25;
  K = 3*sqrt(3)/2;
  
  % Activity phase parameters
  alpha = 0.6;
  beta = 1;
  rho = 0.05;

  try
    % Try and see if there any parameters we should load
    params = r.info.params;
    
    fn = fieldnames(params)
    
    for i = 1:numel(fn)
      cmd = sprintf('%s = %d;', fn{i}, getfield(params,fn{i}));
      fprintf('Setting: %s = %d\n', fn{i}, getfield(params,fn{i})); 
      eval(cmd);
    end
  catch e
    getReport(e)
    keyboard
  end
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % In Grimbert and Cang the positions in SC represents arbors, and
  % they can move around.
  %
  % What I do to repliate their results:
  % Place a dense grid on the SC
  % Place the arbors on this grid
  % Move all arbor points coordinates to r.SCap and r.SCml, remove grid
  % Update connection matrix
  %
  
  if(~exist('useAltP') | isempty(useAltP)) 
    useAltP = true;
  end

  if(~exist('Pfactor'))
    Pfactor = 1;
  end
  
  if(useAltP)
    fprintf('Using altP: P = 1/(f+r).^%g',Pfactor)
  else
    disp('Using originalP: P = 1/(f*r)')
  end
  
  
  if(~exist('ignoreGradients'))
    ignoreGradients = false; % Debug mode, to use Grimbert & Cang
                             % default WT P-functions (regardless of phenotype)
  end

  debugGradients = false;
  
  if(debugGradients)
    makeTestGradients()
    return
  end
    
  r.info.useAltP = useAltP;
  r.info.Pfactor = Pfactor;      
  r.info.ignoreGradients = ignoreGradients;
  
  
  % Place SC neurons as a grid. There is no real notion of SC cells in the
  % model, since the arbor points move around in the SC.
  nGridSC = 100;
  r.nSC = nGridSC^2;
  r.placeSCregularGrid();
  r.loadGradients(r.phenotype);

  r.nSteps = 500;
  placeArbors();

  stepGC();
  % r.plotFigures = 1;
  
  % Some book-keeping necessary
  r.totalWeightRGC = nArbors*ones(r.nRGC,1);
  r.totalWeightSC = ones(r.nSC,1);
  r.saveState();
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % 1. Convert the gradients to probabilities, along the NT axis

  function placeArbors()
    
    if(useAltP)
      
      if(allowRenameSimName)
        r.simName = sprintf('%s-altP',r.simName);
      end
      
      disp('Using alternative P formulation')
      % Alternative formulation for P
      P = 1./(kron(r.RGCEphA,transpose(r.SCephrinA)) ...
              + kron(r.RGCephrinA,transpose(r.SCEphA)));

      if(Pfactor ~= 1)
        P = P.^Pfactor;
        
        if(allowRenameSimName)
          r.simName = sprintf('%s-altP-%g',r.simName, Pfactor);
        end
      end
      
    else
      if(allowRenameSimName)
        r.simName = sprintf('%s-origP',r.simName);
      end
      
      disp('Using original P formulation')
      % Original formulation for P
      P = 1./(kron(r.RGCEphA,transpose(r.SCephrinA)) ...
              .* kron(r.RGCephrinA,transpose(r.SCEphA)));
      
    end

    if(all(P == inf))
      disp('EphA or ephrinA are zero, giving infinite probabilities')
      disp('Forcing P = 1 for all cases')
      P(:) = 1;
    end
    
    % Normalise so max P is 1
    % P = P / max(P(:));
    
    % Note RGC are on row, and SC on column
    % So we need to normalise the rows to 1
    P = P ./ repmat(max(P,[],2),1,size(P,2));
    
    % This overwrites the gradient P function with Grimbert and
    % Cangs WT P value
    if(ignoreGradients)
      disp('WARNING: You are using Grimbert and Cang P-values')
      beep
      r.simName = 'GrimbertCang-original-P-WT';
      bFT = 0.35;
      bFN = 1.40;
      bRT = -0.15;
      bRN = 0.45;
      rF = 20;
      rR = -20; % Should be negative??
      
      try
      BF = (bFT - bFN)*kron(r.RGCnt,ones(1,r.nSC)) + bFN - kron(ones(r.nRGC,1),transpose(r.SCap));
      BR = (bRT - bRN)*kron(r.RGCnt,ones(1,r.nSC)) + bRN - kron(ones(r.nRGC,1),transpose(r.SCap));
      catch e
        getReport(e)
        keyboard
      end
      
      % !!! Work in progress, trying to figure out what rF and rR
      % is set to in their model, not obvious
      PF = 1./(1 + exp(-rF*BF));
      PR = 1./(1 + exp(-rR*BR));
      
      P = PF .* PR;
      
      % Normalise ? (Does not make a big difference, peak is at
      % around 0.99 for GC)
      P = P ./ repmat(max(P,[],2),1,size(P,2));
      
    end
    
    if(0)
      % Make figure for debugging
      figure
      [xP,yP] = meshgrid(linspace(0,1,20),linspace(0,1,20));
      zP = griddata(kron(r.RGCnt,ones(1,r.nSC)), ...
                    kron(ones(r.nRGC,1),transpose(r.SCap)), ...
                    P,xP,yP);
      surf(xP,yP,zP);
      xlabel('Nasal - Temporal')
      ylabel('Anterior - Posterior')
      zlabel('P (unscaled)')
    end
    % P = P/sum(P(:)); % Normalisation is done later

    gridML = unique(r.SCml);

    % We need to establish the borders of the ML, since r.SCml has
    % been rescaled, and we need the unscaled borders to place new
    % arbors.
    [minML,maxML,minAP,maxAP] = SCrange();

    arborML = zeros(r.nRGC,nArbors);
    arborAP = zeros(r.nRGC,nArbors);
  
    % 2. For each RGC neuron, place arbors
    for i = 1:r.nRGC
  
      % Randomize the ML position of the axon
      % Second line makes sure the axon is within the SC
      axonML = (1-r.RGCdv(i)) + sigmaE*randn(1);
      axonML = max(min(axonML,1), 0);
      
      r.RGCml(i) = axonML;
      
      idxCandidate = zeros(nArborsInit,1);
      
      for iArbor = 1:nArborsInit
        picked = false;
        
        trials = 0;
        
        while(~picked)
          trials = trials + 1;
          
          thisArborML = axonML + sigmae*randn(1);
        
          [~,MLidx] = min(abs(gridML-thisArborML));
          SCidx = find(r.SCml == gridML(MLidx));
        
          %P = 1./(r.RGCEphA(i)*r.SCephrinA(SCidx) ...
          %        + r.RGCephrinA(i)*r.SCEphA(SCidx));
        
          randIdx = ceil(rand(1)*numel(SCidx));
          Palt = P(i,SCidx(randIdx));

          if(rand(1) < Palt)
            picked = true;
          end
          
          % Check that arbor inside SC mask
          % The scaling with maxML and minML is since we got some
          % padding around the mask image. This is usually taken
          % care of by rescaling the coordinates later, but here we
          % cant do that.
          AP = r.SCap(SCidx(randIdx));
          if(picked & isempty(r.filterSCPosition(minAP+AP*(maxML-minML),...
                                                 minML+thisArborML*(maxML-minML))))
            % Outside of SC, reject
            picked = false;

            % fprintf('Rejected: %d %d\n', AP, thisArborML)
          end
          
          if(trials > 1e4)
            disp('runGrimbertCang: Many many trials, something is seriously wrong.')
            keyboard
          end
        end
        
        idxCandidate(iArbor) = SCidx(randIdx);
      end
      
      % Decide which of the seven potential arbors we keep
      % P = 1./(r.RGCEphA(i)*r.SCephrinA(idxCandidate) ...
      %        + r.RGCephrinA(i)*r.SCEphA(idxCandidate));
      
      Parbor = transpose(P(i,idxCandidate));
      
      % Add a bit of noise, as they do in their paper
      Parbor = Parbor .* (1 + sigmaNs*randn(nArborsInit,1));
      Parbor = Parbor / sum(Parbor);
      
      [~,idxKeep] = sort(Parbor,'descend');
      SCidxKeep = idxCandidate(idxKeep(1:nArbors));
      
      arborML(i,:) = r.SCml(SCidxKeep);
      arborAP(i,:) = r.SCap(SCidxKeep);  
      
      % keyboard
    end
    
    % Now that we got the points, update the data structures

    % Clear the Eph and ephrin variables, since the model moves the
    % arbors around, so these will contain information not matching
    % the new locations.
    r.RGCEphA = [];
    r.RGCephrinA = [];
    r.SCephrinA = [];
    r.SCEphA = [];
    
    r.RGCEphB = [];
    r.RGCephrinB = [];
    r.SCephrinB = [];
    r.SCEphB = [];

      
    r.SCap = zeros(nArbors*r.nRGC,1);
    r.SCml = zeros(nArbors*r.nRGC,1);
    r.nSC = nArbors*r.nRGC;
    
    r.presynapticConnections = zeros(1,r.nSC);
    r.presynapticWeight = ones(1,r.nSC);
    r.numPresynapticConnections = ones(1,r.nSC);
    
    idxSC = 1;
    for i = 1:r.nRGC
      for j = 1:nArbors
        r.presynapticConnections(1,idxSC) = i;
        r.SCap(idxSC) = arborAP(i,j);
        r.SCml(idxSC) = arborML(i,j);
        
        idxSC = idxSC + 1;
      end
    end
    
  end
  
  
  % End of chemoaffinity phase

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Dynamic (activity) phase begins

  function stepGC(nSteps)

    if(~exist('nSteps'))
      nSteps = r.nSteps - r.curStep;
    end
    
    % Since RGC do not move, precompute all W

    % In the original G&C code they are doing something seriously weird
    % normalisation, so that the smallest value is 1.
   
    NT = kron(r.RGCnt,[1;1;1]);
    DV = kron(r.RGCdv,[1;1;1]);
    
    % Connections are set up so that NT(1),DV(1) proj to AP(1),ML(1)
    if(r.plotFigures)
      figure
      plot(r.SCap,r.SCml,'k.')
    end
      
    xAct = zeros(size(r.SCap));
    yAct = zeros(size(r.SCml));
    
    for iStep = 1:nSteps

      fprintf('Iteration %d/%d\n', iStep, r.nSteps)
      
      % Calculate the activity contribution

      for i = 1:r.nSC
        W = exp(-((NT-NT(i)).^2 + (DV-DV(i)).^2)/rho^2);
        [Fx,Fy] = Ffunc(r.SCap-r.SCap(i),r.SCml-r.SCml(i));
        xAct(i) = alpha*sum(Fx.*W);
        yAct(i) = alpha*sum(Fy.*W);
      end
      
      % Calculate the density contribution
      [gradientX,gradientY,density] = estimateDensityGradient();

      % Update arbor positions in SC
      r.SCap = r.SCap + dt*(alpha*xAct - beta*gradientX);
      r.SCml = r.SCml + dt*(alpha*yAct - beta*gradientY);

      % Prevent exiting 0-1 square
      r.SCap = max(0,min(1,r.SCap));
      r.SCml = max(0,min(1,r.SCml));
      
      if(nnz(r.SCap < 0 | r.SCap > 1 | r.SCml < 0 | r.SCap > 1))
        disp('Outside range!')
        keyboard
      end
      
      if(0)
        image(linspace(0,1,512),linspace(0,1,512),density*50),hold on
        plot(r.SCap,r.SCml,'k.'), hold off
        drawnow
      end
    end
    
    r.curStep = r.curStep + nSteps;
    
  end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function [Fx,Fy] = Ffunc(dx,dy)
    
    % G&C scale by dr^2 in their regular grid, we have SC area
    % approx 0.56 with the current shape, so our scaling is approximately:
    scaleFactor = 0.56/r.nSC;
    
    d2 = dx.^2+dy.^2;
    F = (K/theta)*(1-d2/(theta^2));
    F(d2 > theta^2) = 0;
    Fx = F.*dx*scaleFactor;
    Fy = F.*dy*scaleFactor;
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
  function [gradientX,gradientY,density] = estimateDensityGradient()
  
    nGrid = 512;
    smoothingRadius = 0.06;
    
    useIdx = nGrid/2+1:nGrid*3/2;    
    
    if(isempty(densityNormalisation))
      [minML,maxML,minAP,maxAP] = SCrange();
      
      % Calculate normalisation matrix for density, based on SC shape
      
      % Need to pick coordinate range for ML that will give valid
      % points, also need to make sure X and Y scales are the same
      edges = linspace(0,1,nGrid);
      
      [Xgrid,Ygrid] = meshgrid(edges,edges);
      de = edges(2)-edges(1);
      
      mask = 1/(2*pi*smoothingRadius^2) ...
             * exp(-((Xgrid-0.5).^2 + (Ygrid-0.5).^2)/(2*smoothingRadius^2));
      
      % maxML and minML is scaling necessary to account for the
      % padding around the image. Normally we rescale the resulting
      % coordinates to the range 0-1, but that does not work here.
      [~,~,~,~,SCmask] = r.filterSCPosition(minAP+Xgrid*(maxML-minML),...
                                            minML+Ygrid*(maxML-minML));
      
      % The SCmaskPad makes sure the neurons stays within the SC,
      % the exp factor is so that we have a gradient outside of the
      % SC, so that any neurons lost there will find their way
      % back. If the entire padding was the same value, there would
      % be no gradient telling the neurons where to go towards.
      SCmaskPad = (1-SCmask).*exp(((Xgrid-0.5).^2 + (Ygrid-0.5).^2))*2;
      SCmask(nGrid*2,nGrid*2) = 0;
      mask(nGrid*2,nGrid*2) = 0;
      
      cte = real(ifft2(fft2(SCmask).*fft2(mask)));
      densityNormalisation = 1./cte(useIdx,useIdx);
    end
      
    % Calculate gradient

    % This makes a 2D histogram of the data
    % Note SCap and SCml must be in range 0 to 1.
    try
      AP = max(1,min(nGrid,ceil(r.SCap*nGrid)));
      ML = max(1,min(nGrid,ceil(r.SCml*nGrid)));
      density = accumarray([ML AP],1,[nGrid nGrid]);
    catch e
      getReport(e)
      keyboard
    end
      
    density(2*nGrid,2*nGrid) = 0;
    density = real(ifft2(fft2(density).*fft2(mask)));
    density = density(useIdx,useIdx).*densityNormalisation;
    density = density / mean(density(:));
    
    % keyboard
    
    % Add padding around
    density = density + SCmaskPad;

    [dGradX,dGradY] = gradient(density); 
    % gradient assumes spacing of 1, so we should actually multiply
    % the gradients by 1/de (spacing of lattice)
       
    gradientX = interp2(Xgrid,Ygrid,dGradX,r.SCap,r.SCml,'linear');
    gradientY = interp2(Xgrid,Ygrid,dGradY,r.SCap,r.SCml,'linear');

  end
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % This is just to visualise how the gradients would look
  % Not normally run.
  
  function makeTestGradients()

    r.nRGC = 50;
    r.nSC = 51;
    
    r.RGCnt = transpose(linspace(0,1,r.nRGC));
    r.RGCdv = 0.5*ones(size(r.RGCnt));
    
    r.SCap = transpose(linspace(0,1,r.nSC));
    r.SCml = 0.5*ones(size(r.SCap));
    
    % r.makeNaiveGradients();
    r.loadGradients(r.phenotype);
    
    PF = (1./kron(r.RGCEphA,transpose(r.SCephrinA))).^1;
    PR = (1./kron(r.RGCephrinA,transpose(r.SCEphA))).^1;
  
    Palt = 1./(kron(r.RGCEphA,transpose(r.SCephrinA)) ...
               + kron(r.RGCephrinA,transpose(r.SCEphA)));

    % Ooops, we want NT on X-axis, AP on Y-axis
    PF = transpose(PF);
    PR = transpose(PR);
    Palt = transpose(Palt);
    
    
    % Normalise
    PR = PR/sum(PR(:));
    PF = PF/sum(PF(:));
  
    %Palt = Palt / sum(Palt(:));
    Palt = Palt ./ repmat(max(Palt,[],1),size(Palt,1),1);
    
    P = PF.*PR;

    %P = P / sum(P(:));
    P = P ./ repmat(max(P,[],1),size(P,1),1);
    
    
    if(0)
      % All in one subplot
      figure
      subplot(2,2,1)
      imagesc(PF)
      title('P_{forward}')
      colorbar
      axis equal
      subplot(2,2,2)
      imagesc(PR)
      colorbar
      title('P_{reverse}')
      axis equal
      subplot(2,2,3)
      imagesc(P)
      colorbar
      title('P_{forward} \cdot P_{reverse}')
      axis equal
      subplot(2,2,4)
      imagesc(Palt)
      title('P_{FR}')
      colorbar
      axis equal
    else
      % Separate plots, better for paper.
      
      bFT = 0.35;
      bFN = 1.40;
      bRT = -0.15;
      bRN = 0.45;
      rF = 20;
      rR = -20; % Should be negative??
      
      try
      BF = (bFT - bFN)*kron(r.RGCnt,ones(1,r.nSC)) + bFN - kron(ones(r.nRGC,1),transpose(r.SCap));
      BR = (bRT - bRN)*kron(r.RGCnt,ones(1,r.nSC)) + bRN - kron(ones(r.nRGC,1),transpose(r.SCap));
      catch e
        getReport(e)
        keyboard
      end
      
      % !!! Work in progress, trying to figure out what rF and rR
      % is set to in their model, not obvious
      PFGC = 1./(1 + exp(-rF*BF));
      PRGC = 1./(1 + exp(-rR*BR));
      
      PGC = PFGC .* PRGC;
          
      % PGC = PGC ./ sum(PGC(:));
      
      % To get axes right way around
      PGC = transpose(PGC);
      
      PGC = PGC ./ repmat(max(PGC,[],1),size(PGC,1),1);      
      

      
      figure
      colormap(flipud(colormap('gray')))
      imagesc([0 1],[0 1],PGC,[0 1])
      colorbar
      set(gca,'xtick',[0.1 0.9],'xticklabel', {'N','T'})
      set(gca,'ytick',[0.1 0.9],'yticklabel', {'A','P'})      
      set(gca,'ticklength',[0 0],'fontsize',24)

      % xlabel('Nasal - Temporal','fontsize',30)
      % ylabel('Anterior - Posterior','fontsize',30)
      
      title('P_{Grimbert}','fontsize',36)
      set(gca,'ydir','normal','fontsize',30)
      axis equal
      axis tight
      %box off
      saveas(gcf,'FIGS/GrimbertCang/PGrim.pdf')
      
      figure
      colormap(flipud(colormap('gray')))
      imagesc([0 1],[0 1],P,[0 1])
      colorbar
      set(gca,'xtick',[0.1 0.9],'xticklabel', {'N','T'})
      set(gca,'ytick',[0.1 0.9],'yticklabel', {'A','P'})      
      set(gca,'ticklength',[0 0],'fontsize',24)
      % xlabel('Nasal - Temporal','fontsize',30)
      % ylabel('Anterior - Posterior','fontsize',30)
      
      title('P_{forward} \cdot P_{reverse}','fontsize',36)
      set(gca,'ydir','normal','fontsize',30)
      axis equal
      axis tight
      %box off
      saveas(gcf,'FIGS/GrimbertCang/PF_PR.pdf')
      
      figure
      colormap(flipud(colormap('gray')))
      imagesc([0 1],[0 1],Palt, [0 1])
      colorbar
      set(gca,'xtick',[0.1 0.9],'xticklabel', {'N','T'})
      set(gca,'ytick',[0.1 0.9],'yticklabel', {'A','P'})      
      set(gca,'ticklength',[0 0],'fontsize',24)     
      % xlabel('Nasal - Temporal','fontsize',30)
      % ylabel('Anterior - Posterior','fontsize',30)
      title('P_{Alt}','fontsize',36)
      set(gca,'ydir','normal','fontsize',30)
      axis equal
      axis tight
      %box off
      saveas(gcf,'FIGS/GrimbertCang/Palt.pdf')      

      
      figure
      colormap(flipud(colormap('gray')))

      Pfactor = 3;
      
      scaleFact = 1;
      
      PFac = 1./(kron(r.RGCEphA,transpose(r.SCephrinA)) ...
               + scaleFact*kron(r.RGCephrinA,transpose(r.SCEphA))).^Pfactor;

      PFac = transpose(PFac);
      % PFac = Palt.^Pfactor;
      
      PFac = PFac ./ repmat(max(PFac,[],1),size(PFac,1),1);            
      % PFac = PFac ./ sum(PFac(:));
      imagesc([0 1],[0 1],PFac,[0 1])
      colorbar
      set(gca,'xtick',[0.1 0.9],'xticklabel', {'N','T'})
      set(gca,'ytick',[0.1 0.9],'yticklabel', {'A','P'})      
      set(gca,'ticklength',[0 0],'fontsize',24)     
      % xlabel('Nasal - Temporal','fontsize',30)
      % ylabel('Anterior - Posterior','fontsize',30)
      title(sprintf('P_{Alt %g}',Pfactor),'fontsize',36)
      set(gca,'ydir','normal','fontsize',30)
      axis equal
      axis tight
      %box off
      saveas(gcf,sprintf('FIGS/GrimbertCang/Palt-Pfactor-%g.pdf',Pfactor))      
      
    end
      
      

  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function [minML,maxML,minAP,maxAP] = SCrange()
    
    x = linspace(0,1,1000);
    [APall,MLall] = meshgrid(x,x);
    
    [APfilt,MLfilt] = r.filterSCPosition(APall,MLall);
    minML = min(MLfilt);
    maxML = max(MLfilt);
    
    minAP = min(APfilt);
    maxAP = max(APfilt);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  end