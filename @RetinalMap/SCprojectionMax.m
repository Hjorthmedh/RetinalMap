function [SCmaxNT,SCmaxDV] = SCprojectionMax(obj,useTableFlag)

  if(~exist('useTableFlag'))
    useTable = true;
  end

  if(useTableFlag)
  
    % Make sure the variable is initialised
    assert(~isempty(obj.presynapticConnections));

    SCmaxNT = zeros(obj.nSC,1);
    SCmaxDV = zeros(obj.nSC,1);
    
    for i = 1:obj.nSC
      nCon = obj.numPresynapticConnections(i);
      [~,idx] = max(obj.presynapticWeight(1:nCon,1));
      
      SCmaxNT(i) = obj.RGCnt(obj.presynapticConnections(idx,i));
      SCmaxDV(i) = obj.RGCdv(obj.presynapticConnections(idx,i));
      
    end
  
  else

    assert(~isempty(obj.connectionMatrix))
    
    % Use connection matrix info
    % One column per RGC, one row per SC
    [~,maxIdx] = max(obj.connectionMatrix,[],2);
    
    SCmaxNT = obj.RGCnt(maxIdx);
    SCmaxDV = obj.RGCdv(maxIdx);
    
  end
    
end