% There is also a complementary function, SCprojectionMAX, that
% does SC -> RGC

function [RGCmaxAP,RGCmaxML] = RGCprojectionMax(obj,useTableFlag)

  if(~exist('useTableFlag'))
    useTable = true;
  end

  if(useTableFlag)
  
    % Make sure the variable is initialised
    assert(~isempty(obj.presynapticConnections));

    obj.convertConnectionTables('pre2post');
  
    RGCmaxAP = zeros(obj.nRGC,1);
    RGCmaxML = zeros(obj.nRGC,1);
  
    for i = 1:obj.nRGC
      nCon = obj.numPostsynapticConnections(i);
      [~,idx] = max(obj.postsynapticWeight(1:nCon,i));
      
      RGCmaxAP(i) = obj.SCap(obj.postsynapticConnections(idx,i));
      RGCmaxML(i) = obj.SCml(obj.postsynapticConnections(idx,i));    
    end
  
  else    

    assert(~isempty(obj.connectionMatrix))
    
    % Use connection matrix info
    % One column per RGC, one row per SC
    [~,maxIdx] = max(obj.connectionMatrix,[],1);
    
    RGCmaxAP = obj.SCap(maxIdx);
    RGCmaxML = obj.SCml(maxIdx);
    
  end  

end