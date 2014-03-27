% Running ten repeats of the model comparison for all phenotypes

for runID = 1:10
  ComparePhenotypeModelling('WT','run',runID,false,'Gierer2D');
  ComparePhenotypeModelling('Isl2homozygous','run',runID,false,'Gierer2D');
  ComparePhenotypeModelling('Isl2heterozygous','run',runID,false,'Gierer2D');
  ComparePhenotypeModelling('TKO','run',runID,false,'Gierer2D');
  ComparePhenotypeModelling('Math5','run',runID,false,'Gierer2D');

  % ComparePhenotypeModelling('ephrinA2mm','run',runID,false,'Gierer2D');
  % ComparePhenotypeModelling('ephrinA5mm','run',runID,false,'Gierer2D');
  % ComparePhenotypeModelling('ephrinA2mmA5mm','run',runID,false,'Gierer2D');

end