% Running ten repeats of the model comparison for all phenotypes

for runID = 1001:1010
  ComparePhenotypeModelling('WT','run',runID,false,'koulakov');
  ComparePhenotypeModelling('Isl2homozygous','run',runID,false,'koulakov');
  ComparePhenotypeModelling('Isl2heterozygous','run',runID,false,'koulakov');
  ComparePhenotypeModelling('TKO','run',runID,false,'koulakov');
  ComparePhenotypeModelling('Math5','run',runID,false,'koulakov');

  % ComparePhenotypeModelling('ephrinA2mm','run',runID,false,'koulakov');
  % ComparePhenotypeModelling('ephrinA5mm','run',runID,false,'koulakov');
  % ComparePhenotypeModelling('ephrinA2mmA5mm','run',runID,false,'koulakov');

end