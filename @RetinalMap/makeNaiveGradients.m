% Example of varargin: 'partialGradientRetinaA',0.3, 'withNoise'
% Note that withNoise should be last, otherwise the noise gets
% removed from the clamped part of the trace


function makeNaiveGradients(obj,varargin)

  disp('Setting up naive gradients')
  
  assert(strcmp(obj.eyeType,'disk'));

  % Set up default gradients first, then allow them to be modified

  % Forward signaling
  % Reber EphA gradient: (use option ReberWT to get this now)
  %  obj.RGCEphA = 2.59/3.64*exp(-2.3*(1-obj.RGCnt)) + 1.05/3.64;
  disp('Not using Reber gradients ---- test.')
  obj.RGCEphA = exp(obj.RGCnt-1);
  obj.RGCEphB = exp(obj.RGCdv-1);

  % Reverse signaling
  obj.RGCephrinA = exp(-obj.RGCnt);
  obj.RGCephrinB = exp(-obj.RGCdv);
  
  % Forward signaling
  obj.SCephrinA = exp(obj.SCap-1);
  obj.SCephrinB = exp(-obj.SCml);
  
  % Reverse signaling
  obj.SCEphA = exp(-obj.SCap);
  obj.SCEphB = exp(obj.SCml-1);
  
  % keyboard
 
  if(nargin > 0)
    ctr = 1;
    
    while(ctr <= nargin - 1) % obj is counted by nargin
  
      switch(varargin{ctr})
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        case 'withNoise'
          
          if(ctr < nargin-1)
            disp(['You should put withNoise last if you want the ' ...
                  'clamped part of the gradients to have noise.'])
          end

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

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        case 'partialGradientRetinaA'
          % This requires two additional arguments, which gradient, 
          % and how far.
      
          if(ctr + 1 > nargin)
            disp('makeNaiveGradients: partialGradientRetinaA fraction')
            keyboard
          end
          
          % 'RGCEphA', 'RGCephrinA'
          clampFraction = varargin{ctr+1};
          
          % Mark that we read an additional arguments
          ctr = ctr + 1;
          
          RGCwidth = max(obj.RGCnt(:))-min(obj.RGCnt(:));
          minNT = min(obj.RGCnt(:));
          
          switch(obj.typeFlag)
            case 1
              % Single gradient, only modify the retinal EphA
                      
              idx = find(obj.RGCnt < minNT + RGCwidth*clampFraction);
              clampVal = max(obj.RGCEphA(idx));
              obj.RGCEphA(idx) = clampVal;

              fprintf('Clamping RGC EphA, fraction %.3f, value %.3f\n', ...
                      clampFraction, clampVal)             
              
            case {2,3,5}
              % Double gradient (or servo), modify retinal EphA and ephrinA

              idx = find(obj.RGCnt < minNT + RGCwidth*clampFraction);
              clampVal = max(obj.RGCEphA(idx));
              obj.RGCEphA(idx) = clampVal;

              fprintf('Clamping RGC EphA, fraction %.3f, value %.3f\n', ...
                      clampFraction, clampVal)              
              
              idx = find(obj.RGCnt > minNT + RGCwidth*(1-clampFraction));
              clampVal = max(obj.RGCephrinA(idx));
              obj.RGCephrinA(idx) = clampVal;

              fprintf('Clamping RGC ephrinA, fraction %.3f, value %.3f\n', ...
                      clampFraction, clampVal)
              
          end
          
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
        
        case 'ReberWT'
          
          % Replace the EphA gradient in the retina with the
          % gradient from Reber et al 2004
          % We rescaled the amplitude of Ra so that the maximum is 1
          % as before (/3.64), this means all Isl2 values are rescaled also
          % At x = 0 we want Ra = 1 for WT, thus the scale factor is
          % 0.26*exp(2.3) + 1.05 = 3.64
        
          
          obj.RGCEphA = 2.59/3.64*exp(-2.3*(1-obj.RGCnt)) + 1.05/3.64;
          
        case 'ReberIsl2heterozygot'
          
          obj.RGCEphA = 2.59/3.64*exp(-2.3*(1-obj.RGCnt)) + 1.05/3.64;          
          % Half of the neurons get an EphA3 boost (heterozygot)

          % heterozygot = (1.98-1.05)/3.64 = 0.2555 
          idx = find(rand(size(obj.RGCEphA)) < 0.5);
          obj.RGCEphA(idx) = obj.RGCEphA(idx) + 0.255;
          
          obj.Isl2PositiveRGC = idx;
          
          
        case 'ReberIsl2homozygot'
          obj.RGCEphA = 2.59/3.64*exp(-2.3*(1-obj.RGCnt)) + 1.05/3.64;        
          % Half of the neurons get an EphA3 boost (homozygot)

          % homozygot = (2.91-1.05)/3.64 = 0.511
          idx = find(rand(size(obj.RGCEphA)) < 0.5);
          obj.RGCEphA(idx) = obj.RGCEphA(idx) + 0.511;

          obj.Isl2PositiveRGC = idx;          
          
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        case 'ephrin-A5-KO'
          
          % If we are using double gradients, zero out the retinal
          % ephrinA5
          obj.RGCephrinA = 0*obj.RGCephrinA;
          
          % Make the gradient shallower, only ephrin-A2 remaining
          obj.SCephrinA = exp(obj.SCap/2-1);
          
        case 'ephrin-A2-KO'
        
          % No ephrin-A2 or A6 in retina
          
          % Remaining ephrin-A5 in SC gives only partial gradient
          
          obj.SCephrinA = max(0,exp(obj.SCap-1)-0.5);
          
        case 'ephrin-A2-A5-KO'

          obj.RGCephrinA = 0*obj.RGCephrinA;
          obj.SCephrinA = 0*obj.SCephrinA;
          
          
        % These are the gradients from the Tsigankov Koulakov
        % 2006 paper
        
        case 'TK2006-WT'

          disp('Tsigankov-Tripplet 2006 WT gradients')
          % Reverse gradient down scaled in their model
          % Note that their model is a typeFlag = 1 model,
          % the gradients are only used for masking
          obj.SCEphA = exp(-obj.SCap) * exp(-1);
          obj.RGCephrinA = exp(-obj.RGCnt) * exp(-1);

          assert(obj.typeFlag == 1);
          
          obj.RGCEphA = max(0,obj.RGCEphA - obj.RGCephrinA);
          obj.SCephrinA = max(0,obj.SCephrinA - obj.SCEphA);
          
        case 'TK2006-ephrinA2KO'

          disp('Tsigankov-Tripplet 2006 ephrinA2 KO gradients')          
          
          SCephrinA5 = 0.66*exp(3*(obj.SCap-1)) +0.15;
          obj.SCephrinA = SCephrinA5;
          
          % Reverse gradient down scaled in their model
          % and only used for masking
          obj.SCEphA = exp(-obj.SCap) * exp(-1);
          obj.RGCephrinA = exp(-obj.RGCnt) * exp(-1);

          assert(obj.typeFlag == 1);          
          
          % Calculate masking
          obj.RGCEphA = max(0,obj.RGCEphA - obj.RGCephrinA);
          obj.SCephrinA = max(0,obj.SCephrinA - obj.SCEphA);
          
        case 'TK2006-ephrinA5KO'
          
          disp('Tsigankov-Tripplet 2006 ephrinA5 KO gradients')          
          
          SCephrinA = exp(obj.SCap-1);
          SCephrinA5 = 0.66*exp(3*(obj.SCap-1)) + 0.15;
          SCephrinA2 = SCephrinA - SCephrinA5;
          obj.SCephrinA = SCephrinA2;

          % Reverse gradient down scaled in their model          
          % and only used for masking
          obj.SCEphA = exp(-obj.SCap) * exp(-1);          
          obj.RGCephrinA = exp(-obj.RGCnt) * exp(-1);
                    
          assert(obj.typeFlag == 1);          
          
          % Calculate masking
          obj.RGCEphA = max(0,obj.RGCEphA - obj.RGCephrinA);
          obj.SCephrinA = max(0,obj.SCephrinA - obj.SCEphA);
        
        case 'TK2006-ephrinA2A5KO'
          
          disp('Tsigankov-Tripplet 2006 ephrinA2A5 KO gradients')          

          % Reverse gradient down scaled in their model
          % and only used for masking
          obj.SCEphA = exp(-obj.SCap) * exp(-1);
          obj.RGCephrinA = exp(-2)*ones(size(obj.RGCephrinA));

          obj.SCephrinA = zeros(size(obj.SCephrinA));
          
          assert(obj.typeFlag == 1);          
          
          % Calculate masking
          obj.RGCEphA = max(0,obj.RGCEphA - obj.RGCephrinA);
          obj.SCephrinA = max(0,obj.SCephrinA - obj.SCEphA);
          
          
        case 'TK2011'
          
          % These are the default Koulakov-Tripplet gradients,
          % unscaled. Tripplet uses TN and VD axis
          X = 1-obj.RGCnt;
          Y = 1-obj.RGCdv;
          
          % His reverse gradients have the same direction in the
          % supplementary material...
          
          obj.RGCEphA = exp(-X);
          obj.RGCEphB = exp(1-Y);
          obj.RGCephrinA = exp(1-X); % WRONG?!
          obj.RGCephrinB = exp(-Y); % WRONG!?
          
          x = 1-obj.SCap;
          y = obj.SCml;
          
          obj.SCephrinA = exp(1-x);
          obj.SCephrinB = exp(1-y);
          obj.SCEphA = exp(-x); 
          obj.SCEphB = exp(-y);
          
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        otherwise
          
          fprintf('Unknown gradient parameter: %s\n', varargin{ctr})
        
      end
      
      ctr = ctr + 1;
    end
  end
    

end