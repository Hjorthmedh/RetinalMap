function initializeRandomGenerator(obj)

  disp('Setting new random seed.')

  if(verLessThan('matlab','7.12'))
    % Edinburgh workaround, support old version of random
    % generator initialisation
    disp('Edinburgh legacy mode')
    seed = sum(1e6*clock);
    s = RandStream.create('mt19937ar', 'seed', seed);
    RandStream.setDefaultStream(s);
  else
    % Default way of generating random seed in Matlab 2012
    rng('default')
    rng('shuffle')
  end

end