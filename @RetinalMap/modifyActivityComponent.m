function modifyActivityComponent(obj,gammaAct,bAct,aAct)

  if(exist('gammaAct') & ~isempty(gammaAct))
    obj.gammaAct = gammaAct;
  end
  
  if(exist('bAct') & ~isempty(bAct))
    obj.bAct = bAct;
  end

  if(exist('aAct') & ~isempty(aAct))
    obj.aAct = aAct;
  end
  
  obj.calculateCU();
  
end