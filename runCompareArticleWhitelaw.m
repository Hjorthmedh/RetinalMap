% Running ten repeats of the model comparison for all phenotypes

for runID = 1001:1010
  ComparePhenotypeModelling('WT','run',runID,false,'WhiteCow');
  ComparePhenotypeModelling('Isl2homozygous','run',runID,false,'WhiteCow');
  ComparePhenotypeModelling('Isl2heterozygous','run',runID,false,'WhiteCow');
  ComparePhenotypeModelling('TKO','run',runID,false,'WhiteCow');
  ComparePhenotypeModelling('Math5','run',runID,false,'WhiteCow');

  % ComparePhenotypeModelling('ephrinA2mm','run',runID,false,'WhiteCow');
  % ComparePhenotypeModelling('ephrinA5mm','run',runID,false,'WhiteCow');
  % ComparePhenotypeModelling('ephrinA2mmA5mm','run',runID,false,'WhiteCow');

end