function fig = plotQuantificationPanel(obj, report,fontsize)

  if(obj.plotFigures == 0)
    disp('plotParameters: plotFigures = 0, hiding figures!')
    visFlag = 'off';
  else
    visFlag = 'on';
  end
  
  fig = figure('visible',visFlag);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  try
  
    textStr{1} = 'Quantification of map';
    textStr{end+1} = sprintf('Segregation AP: inflection %.3f slope %.3f, Segregation ML: inflection %.3f slope %.3f', ...
                             report.kSegregationAP(1), report.kSegregationAP(2), ...
                             report.kSegregationML(1), report.kSegregationML(2));
  
    textStr{end+1} = sprintf('AP collapse point (0-1): %.2f', ...
                             report.APcollapsePoint);
  
    try
      textStr{end+1} = sprintf(['dFuzz NT = %.3f, kFuzz AP = %.3f, cFuzz95 AP = %.3f\n'...
                                'dFuzz DV = %.3f, kFuzz ML %.3f, cFuzz95 ML = %.3f'],...
                               report.dFuzzNT, report.kFuzzAP, report.cFuzz95AP, ...
                               report.dFuzzDV, report.kFuzzML, report.cFuzz95ML);
    catch
      % No dFuzz measures, skipping.
    end
      
    textStr{end+1} = sprintf('SC covered by 95%% of synapses: %.3f', ...
                             report.SC95coverage);
  
    textStr{end+1} = sprintf('Nasal ectopic: %s, Temporal ectopic: %s', ...
                             yesno(report.nasalEctopic), ...
                             yesno(report.temporalEctopic));

    textStr{end+1} = sprintf(['Retina covered by 5/25/50/75/95 %%' ...
                              'of synapses: %.1f/%.1f/%.1f/%.1f/%.1f %%'], ...
                             100*report.areaContour);
    
    if(~isempty(report.DWgrid))
      % Display DW analysis measures
      
      textStr{end+1} = sprintf('Number of ectopics: %d DW', ...
                               report.DWgrid.num_ectopics);  

      try
        textStr{end+1} = sprintf('Ectopic distance: %.3f +/- %.3f (%.3f +/- %.3f mm) DW', ...
                                 report.DWgrid.ect_dist_mean, ...
                                 report.DWgrid.ect_dist_std, ...
                                 report.DWgrid.ect_dist_mean * obj.SCscale, ...
                                 report.DWgrid.ect_dist_std * obj.SCscale);  
      catch e
        getReport(e)
        % The DW analysis might have crashed, and these variables
        % might not exist. So keep going...
      end
        
      try
        textStr{end+1} = sprintf('Number of crossings: %d DW FTOC (%d)', ...
                                 report.DWgrid.FTOC.num_crossings, ...
                                 report.grid.SCnumCrossings);
  
        textStr{end+1} = sprintf('Number of nodes in subgraph: %d DW FTOC (%d)', ...
                                 report.DWgrid.FTOC.num_nodes_in_subgraph, ...
                                 report.grid.SCnodesLeft);
      catch
        getReport(e)
        % Continue...
      end
  
      try
        textStr{end+1} = sprintf('Projection spread. AP %.3f (%.3fmm), ML %.3f (%.3fmm) DW FTOC (major,minor?)', ...
                                 report.DWgrid.FTOC.overall_dispersion_yrad, ...
                                 report.DWgrid.FTOC.overall_dispersion_yrad * obj.SCscale, ...                               
                                 report.DWgrid.FTOC.overall_dispersion_xrad, ...
                                 report.DWgrid.FTOC.overall_dispersion_xrad * obj.SCscale);
      catch
        disp('Missing lattice dispersion measures.')
      end
    else

      % We didnt run DW analysis, so only display my own measures
      
      textStr{end+1} = sprintf('Number of crossings: (%d)', ...
                               report.grid.SCnumCrossings);
      

      textStr{end+1} = sprintf('Number of nodes in subgraph: (%d)', ...
                               report.grid.SCnodesLeft);
      
    end
      
    textStr{end+1} = sprintf('Major axis angle: %.1f degrees, Major/minor ratio %.2f',...
                             report.superposedProjMajorAng*180/pi, ...
                             report.superposedProjMajorMinorRatio);
  
  catch e
    getReport(e)
    keyboard
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  t = text(0.1,0.5,textStr);
  
  if(exist('fontsize'))
    set(t,'fontsize',fontsize)
  end
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function answer = yesno(flag)

    if(flag)
      answer = 'Yes';
    else
      answer = 'No';
    end
  
  end
  
  axis off
end