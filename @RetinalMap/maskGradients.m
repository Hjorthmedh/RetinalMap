% This function calculates the effective available gradients,
% assuming there is simple masking due to cis interaction.
%
% [R] + [L] <-> [RL], here we assume k = 1 for the reaction
%
% TODO: See DS concern raised March 29, 2012. However, the cis
% interaction takes place over a larger area than the trans
% interaction, so maybe we are justified to do that equilibrium
% first. NEED TO VERIFY THIS.
%
% Tsigankov and Koulakov 2006 implement masking as subtraction, however
% Reber et al 2004, supplementary discussion, argues that if there
% is masking then it is not complete in the nasal retina. They
% motivate it by looking at the effect of EphA4 KO in Isl2-EphA3,
% where if all EphA4 would be inactivated nasally, the KO would have no
% effect nasally.
%
%

function maskGradients(obj, plotFlag)
 
  switch(obj.kMask)
    case 0
      return; % No masking
    case 1
      % Reaction scheme masking, this function, below
    case 2
      obj.maskGradientsSubtractive();
      return
    case 3
      obj.maskGradientsMultiplicative();
      return
  end
  
  disp('!!! Verify that maskGradients equations are correct')
  
  assert(0 <= obj.kMask & obj.kMask <= 1);

  oldRGCEphA = obj.RGCEphA;
  oldRGCephrinA = obj.RGCephrinA;
  oldRGCEphB = obj.RGCEphB;
  oldRGCephrinB = obj.RGCephrinB;
  
  oldSCephrinA = obj.SCephrinA;
  oldSCephrinB = obj.SCephrinB;
  oldSCEphA = obj.SCEphA;
  oldSCEphB = obj.SCEphB;
  
  fprintf('Masking gradients, using kMask = %f\n', obj.kMask)

  % The min below prevents gradients going negative
  boundRGCA = obj.kMask*obj.RGCEphA.*obj.RGCephrinA;
  boundRGCA = min(boundRGCA,min(obj.RGCEphA,obj.RGCephrinA));
  obj.RGCEphA = obj.RGCEphA - boundRGCA;
  obj.RGCephrinA = obj.RGCephrinA - boundRGCA;

  boundRGCB = obj.kMask*obj.RGCEphB.*obj.RGCephrinB;
  boundRGCB = min(boundRGCB,min(obj.RGCEphB,obj.RGCephrinB));
  obj.RGCEphB = obj.RGCEphB - boundRGCB;
  obj.RGCephrinB = obj.RGCephrinB - boundRGCB;
  
  boundSCA = obj.kMask*obj.SCephrinA.*obj.SCEphA;
  boundSCA = min(boundSCA,min(obj.SCephrinA,obj.SCEphA));
  obj.SCephrinA = obj.SCephrinA - boundSCA;
  obj.SCEphA = obj.SCEphA - boundSCA;
  
  boundSCB = obj.kMask*obj.SCephrinB.*obj.SCEphB;
  boundSCB = min(boundSCB,min(obj.SCephrinB,obj.SCEphB));
  obj.SCephrinB = obj.SCephrinB - boundSCB;
  obj.SCEphB = obj.SCEphB - boundSCB;
  
  % Just make sure there are no negative receptor/ligand concentractions
  
  % Stop if any are negative
  try
    assert(nnz(obj.RGCEphA < 0) == 0);
    assert(nnz(obj.RGCEphB < 0) == 0);  
    assert(nnz(obj.RGCephrinA < 0) == 0);  
    assert(nnz(obj.RGCephrinB < 0) == 0);    
      
    assert(nnz(obj.SCephrinA < 0) == 0);  
    assert(nnz(obj.SCephrinB < 0) == 0);    
    assert(nnz(obj.SCEphA < 0) == 0);    
    assert(nnz(obj.SCEphB < 0) == 0);      
  catch e
    disp('Masking resulted in negative gradients!')
    getReport(e)
    keyboard
  end    
    
  % Plot
  if(plotFlag | ~obj.noFigs)

    mSize = 15;

    figure
    % Retina A system
    subplot(2,2,1)
    p(1) = plot(obj.RGCnt, oldRGCEphA, '.', ...
                'markersize', mSize, 'color', [1 0.5 0.5]);
    hold on
    p(2) = plot(obj.RGCnt, oldRGCephrinA, '.', ...
                'markersize', mSize, 'color', [0.5 0.5 0.5]);

    p(3:4) = plot(obj.RGCnt, obj.RGCEphA, 'r.', ...
                  obj.RGCnt, obj.RGCephrinA, 'k.', ...
                  'markersize', mSize);
    hold off
    xlabel('Nasal - Temporal')
    ylabel('Concentration')
    legend(p,'EphA (orig)', 'ephrinA (orig)', ...
             'EphA (masked)', 'ephrinA (masked)')
    title('Retina')
    
    % Retina B system
    subplot(2,2,3)
    p(1) = plot(obj.RGCdv, oldRGCEphB, '.', ...
                'markersize', mSize, 'color', [0.5 0.5 1]);
    hold on
    p(2) = plot(obj.RGCdv, oldRGCephrinB, '.', ...
                'markersize', mSize, 'color', [0.5 0.5 0.5]);

    p(3:4) = plot(obj.RGCdv, obj.RGCEphB, 'b.', ...
                  obj.RGCdv, obj.RGCephrinB, 'k.', ...
                  'markersize', mSize);
    hold off
    xlabel('Dorsal - Ventral')
    ylabel('Concentration')
    legend(p,'EphB (orig)', 'ephrinB (orig)', ...
             'EphB (masked)', 'ephrinB (masked)')

    
    % SC A system
    subplot(2,2,2)
    p(1) = plot(obj.SCap, oldSCephrinA, '.', ...
                'markersize', mSize, 'color', [1 0.5 0.5]);
    hold on
    p(2) = plot(obj.SCap, oldSCEphA, '.', ...
                'markersize', mSize, 'color', [0.5 0.5 0.5]);
    p(3:4) = plot(obj.SCap,obj.SCephrinA, 'r.', ...
                  obj.SCap,obj.SCEphA,'k.', ...
                  'markersize', mSize);
    hold off
    
    xlabel('Anterior - Posterior')
    ylabel('Concentration')
    legend(p, 'ephrinA (orig)', 'EphA (orig)', ...
              'ephrinA (masked)', 'EphA (masked)')

    title('Superior colliculus')

    % SC B system
    subplot(2,2,4)
    p(1) = plot(obj.SCml, oldSCephrinB, '.', ...
                'markersize', mSize, 'color', [0.5 0.5 1]);
    hold on
    p(2) = plot(obj.SCml, oldSCEphB, '.', ...
                'markersize', mSize, 'color', [0.5 0.5 0.5]);
    p(3:4) = plot(obj.SCml,obj.SCephrinB, 'r.', ...
                  obj.SCml,obj.SCEphB,'k.', ...
                  'markersize', mSize);
    hold off
    
    xlabel('Medial - Lateral')
    ylabel('Concentration')
    legend(p, 'ephrinA (orig)', 'EphA (orig)', ...
              'ephrinA (masked)', 'EphA (masked)')
    
    fName = sprintf('%s/%s-maskedGradients.pdf', ...
                    obj.figurePath, obj.simName);
    saveas(gcf,fName,'pdf')
    
  end
  
end