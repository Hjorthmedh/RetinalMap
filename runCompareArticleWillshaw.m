% Running ten repeats of the model comparison for all phenotypes

for runID = [4,6:12,100,101]
  ComparePhenotypeModelling('WT','run',runID,false,0,'Markerinduction');
  ComparePhenotypeModelling('Isl2homozygous','run',runID,false,0,'Markerinduction');
  ComparePhenotypeModelling('Isl2heterozygous','run',runID,false,0,'Markerinduction');
  ComparePhenotypeModelling('TKO','run',runID,false,0,'Markerinduction');
  ComparePhenotypeModelling('Math5','run',runID,false,0,'Markerinduction');

  % ComparePhenotypeModelling('ephrinA2mm','run',runID,false,0,'Markerinduction');
  % ComparePhenotypeModelling('ephrinA5mm','run',runID,false,0,'Markerinduction');
  % ComparePhenotypeModelling('ephrinA2mmA5mm','run',runID,false,0,'Markerinduction');

end