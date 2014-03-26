% What is the minimum fraction of the SC is covered by 95% of all synapses

function fracCovered = calculateSCcoverage(obj,synapseFrac)

  if(~exist('synapseFrac'))
    synapseFrac = 0.95;
  end

  % We start eliminating the SC cells with the fewest synapses
  % How many of those can we take away before we have removed 
  % nSynExclude synapses?
  
  nSynExclude = (1-synapseFrac)*sum(obj.totalWeightSC);
  nSyn = sort(obj.totalWeightSC);

  lastIdxExclude = find(cumsum(nSyn) <= nSynExclude,1,'last');
  
  fracCovered = (obj.nSC-lastIdxExclude)/obj.nSC;

end