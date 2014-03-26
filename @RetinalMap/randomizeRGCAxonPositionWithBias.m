function randomizeRGCAxonPositionWithBias(obj)

  disp('Placing initial axons with initial bias')

  mlRange = [min(obj.SCml) max(obj.SCml)];
  dvRange = [max(obj.RGCdv) min(obj.RGCdv)];

  switch(obj.eyeType)
    case 'sphere'
      assert(false);
      %perfectAP = interp1(dvRange,mlRange, ??? );
    case 'disk'
      perfectAP = interp1(dvRange,mlRange,obj.RGCdv);
  end
      
  obj.RGCml = zeros(obj.nRGC,1);
  
  for i = 1:numel(perfectAP)
    mlPos = perfectAP(i) + obj.RGCmlSpread*randn(1);
    
    % Make sure our random position is within the SC
    while(mlPos < mlRange(1) | mlRange(2) < mlPos )
      mlPos = perfectAP(i) + obj.RGCmlSpread*randn(1);    
    end

    obj.RGCml(i) = mlPos;
    
  end
  
  if(obj.plotFigures)
    % Just a verification plot
    % Compare with Plas, Lopez, Crair 2005
    
    % r1 = [0.1 0.2];
    % r2 = [0.8 0.9];

    r1 = [0 0.1];
    r2 = [0.9 1];
    
    idxD = find(r1(1) < obj.RGCdv & obj.RGCdv < r1(2));
    idxV = find(r2(1) < obj.RGCdv & obj.RGCdv < r2(2));
    idxN = find(r1(1) < obj.RGCnt & obj.RGCnt < r1(2));
    idxT = find(r2(1) < obj.RGCnt & obj.RGCnt < r2(2));
    
    idxAll = {idxD, idxT, idxV, idxN };
    legText = { 'Dorsal injection', ...
                'Temporal injection', ...
                'Ventral injection', ...
                'Nasal injection', ...
                };
    
    col = 'rbcg';
    
    edges = linspace(0,1,50);
    
    pD = load('initialBias/PlasLopezCrair2005-figure4E-dorsal-injection.csv');
    pV = load('initialBias/PlasLopezCrair2005-figure4E-ventral-injection.csv');
    pN = load('initialBias/PlasLopezCrair2005-figure4E-nasal-injection.csv'); 
    pT = load('initialBias/PlasLopezCrair2005-figure4E-temporal-injection.csv');
       
    plas = {pD,pT,pV,pN};
    
    maxY = 0;
    
    figure
    for i = 1:numel(idxAll)
      sp(i) = subplot(2,2,i);
      n = histc(obj.RGCml(idxAll{i}),edges);
      nNorm = n/(sum(n)*(edges(2)-edges(1)));
      p(i) = stairs(edges,nNorm,'color',col(i));
      hold on
      tn = plas{i}(:,1);
      fluor = plas{i}(:,2);
      % We need to normalise it
      w = [0; diff(tn)/2] + [diff(tn/2); 0];
      fluor = fluor / sum(fluor.*w);
      nt = 1-tn;
      plot(nt,fluor,'--','color',col(i))
      title(legText{i})
      xlabel('Medial-Lateral')
      ylabel('Normalised count')
    
      maxY = max([maxY,max(fluor),max(nNorm)]);
      
    end
       
    linkaxes(sp,'xy');
    axis([0 1 0 maxY]);
    
    % legend(p,legText,'location','northwest');

    % title('RGC axon ingrowth bias')
    
    fName = sprintf('%s/%s-RGC-axon-ingrowth-bias-%.2f.pdf', ...
                    obj.figurePath, obj.simName, obj.RGCmlSpread);
    saveas(gcf,fName,'pdf')
    
  end

end