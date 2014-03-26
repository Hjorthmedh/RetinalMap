% Helper function, so that the experimental data from Dans
% segregation fit is not spread throughout the code

function [WT,B2KO,SCsize] = getDanSegregationFits(obj)

  % Table 2-4 in Dan Lyngholms thesis, SC size (age,AP,ML)
  SCsize = [0, 1910, 1800;
            2, 1920, 1960;
            4, 2000, 1900;
            6, 1980, 1850;
            8, 1890, 1840;
            12, 1920, 1890;
            16, 1970, 1980;
            22, 1970, 1960;
            60, 2072, 2008];  
  
  % Table 4-2 in Dan Lyngholms thesis, WT segregation fits
  
  WT.age = [0 2 4 6 8 12 16 22 60];

  WT.kAP = [2749   3.6; 
            2711   1.432;
            753.6  1.331;
            204.9  0.9167;
            212.7  1.086;
            136.4  2.88;
            127.8  2.786;
            112.3  2.195;
            74.88  5.169];
  
  WT.kML = [551.3  2.646;
            661.1  4.72;
            464.2  4.937;
            NaN    NaN;
            231.2  1.374;
            166.2  1.716;
            NaN    NaN;
            119.4  1.995
            NaN    NaN];

  WT.kAPscaled = scaleK(WT.kAP,WT.age,'AP');
  WT.kMLscaled = scaleK(WT.kML,WT.age,'ML');  
  
  % B2KO, Dan Lyngholm thesis, table 4-3

  B2KO.age = [0 2 4 8 12 22 60];
  
  B2KO.kAP = [2354   4.52;
              2934   5.123;
              2464   23.44;
              869.6  1.949;
              406.8  3.164;
              447    3.775;
              404    5.472];

  B2KO.kAPscaled = scaleK(B2KO.kAP,B2KO.age,'AP');
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function dScaled = scaleSCdist(dRaw,age,axisDir)
    % We want the largest axis to be between 0 and 1, so scale by AP size 
          
    SCage = SCsize(:,1);
    SCAP = SCsize(:,2);
    SCML = SCsize(:,3);
    
    switch(axisDir)
      case 'AP'
        scaleFactor = transpose(interp1(SCage,SCAP,age));
      case 'ML'
        scaleFactor = transpose(interp1(SCage,SCML,age));
      otherwise
        fprintf('Unknown axis dir: %s (use ''AP'' or ''ML'')\n', axisDir)
        keyboard
    end

    dScaled = dRaw ./ scaleFactor;
    
  end
   
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function kScale = scaleK(k,age,axisDir)
    
    assert(size(k,1) == numel(age));

    kScale = k;   
    kScale(:,1) = scaleSCdist(k(:,1),age,axisDir);
  
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
end