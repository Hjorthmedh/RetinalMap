function signalNoise = addPoissonNoise(obj,signal,N)

  signalNoise = poissrnd(signal*N)/N;

end