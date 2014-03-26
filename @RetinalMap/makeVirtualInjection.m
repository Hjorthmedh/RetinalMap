function [RGCidx,SCidx] = makeVirtualInjection(obj,SCap,SCml,injRadius)

  if(~exist('injRadius') | isempty(injRadius))
    % !!! Update, from Dans thesis, chap 3.3.1 -- 2.8% of SC, diameter
    % For WT adult, this should give 8.5% of retinal diameter
  
    SCwidth = max(obj.SCap(:))-min(obj.SCap(:));  
    injRadius = SCwidth*0.028/2;
  end
    
  SCidx = find((obj.SCap-SCap).^2 + (obj.SCml-SCml).^2 < injRadius^2);
  RGCidx = obj.findPresynapticRGC(SCidx);
  
end