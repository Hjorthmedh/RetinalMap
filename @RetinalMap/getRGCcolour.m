function RGCcolour = getRGCcolour(obj)

  assert(strcmp(obj.eyeType,'disk'))

  if(ismember('RGCcolour',fieldnames(obj)) & ~isempty(obj.RGCcolour))
    RGCcolour = obj.RGCcolour;
  else
    
    RGCb = (obj.RGCnt - min(obj.RGCnt)) / (max(obj.RGCnt) - min(obj.RGCnt));
    RGCg = (obj.RGCdv - min(obj.RGCdv)) / (max(obj.RGCdv) - min(obj.RGCdv));
  
    RGCcolour = [zeros(size(RGCb)) RGCg RGCb];

  end
    
end