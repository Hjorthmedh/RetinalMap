% This file gets the raw data from Cangs intrinsic imaging
% The data is stored in the map-analysis directory

clear all, close all
% TKO homozygous: 54, 55, 56, 58
% WT: 6, 10, 73, 80
% Isl2 ki: 999 (zebrafinch), 115

dataSet = [58 6 10 73 80 54 55 56 115];

addpath('../map-analysis/lib')
addpath('../map-analysis/cang')

for i = 1:numel(dataSet)

  params = Dgetparams(dataSet(i));

  % load_data uses relative paths, need to be in the cang directory
  % for the function to work...
  curDir = pwd;
  cd('../map-analysis/cang')
  params = load_data(params);
  cd(curDir);
  
  % Piggyback on lattice analysis code
  
  params = find_active_pixels(params);
  params = make_list_of_points(params); 
  
  % HIGH SCATTER POINTS REMOVED
  params = Dremove_high_scatter(params);
  params = make_list_of_points(params); 
  params = Dfind_map_quality(params)

  
  % The data seems to have DV then NT, and ML then AP. For Willshaw
  % grid we send it in as NT,DV and AP,ML... this could be
  % potentially confusing!
  
  NT = params.full_field(:,2);
  NT = (max(NT) - NT)/(max(NT)-min(NT));
  
  DV = params.full_field(:,1);
  DV = (max(DV) - DV)/(max(DV)-min(DV));
  
  AP = params.full_coll(:,2);
  AP = (AP - min(AP))/ (max(AP) - min(AP));
  
  ML = params.full_coll(:,1);
  ML = (ML - min(ML)) / (max(ML) - min(ML));
  
  
  
  figure
  subplot(1,2,1)
  
  % We only want to plot central third (along DV)
  idx = find(1/3 <= DV & DV <= 2/3);
  
  plot(NT(idx),AP(idx),'r.')
  set(gca,'xtick',[0.1 0.9],'xticklabel',{'N','T'})
  set(gca,'ytick',[0.1 0.9],'yticklabel',{'A','P'})
  set(gca,'fontsize',20)
  box off
  
  axis equal
  %  axis tight
  axis([0 1 0 1])


  title(sprintf('ID: %d (%s)', params.id, params.datalabel))
  
  subplot(1,2,2)
  idx = find(1/3 <= NT & NT <= 2/3);
  plot(DV(idx),ML(idx),'r.')
  set(gca,'xtick',[0.1 0.9],'xticklabel',{'D','V'})
  set(gca,'ytick',[0.1 0.9],'yticklabel',{'M','L'})
  set(gca,'fontsize',20)
  box off
  axis equal
  %  axis tight
  axis([0 1 0 1])  


  
  
  fName = sprintf('FIGS/DW-proj-id-%d.pdf', params.id);
  saveas(gcf,fName,'pdf')
  
end

rmpath('../map-analysis/lib')
rmpath('../map-analysis/cang')
  