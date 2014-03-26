% Filter out the points outside our defined region

function [NT,DV,NTpad,DVpad,mask] = filterRetinalPositionsDisk(obj,NT,DV)

  NTc = 0.5; DVc = 0.5; r = 0.5;

  okPos = sqrt((NT-NTc).^2+(DV-DVc).^2) < r;

  idxPad = find(~okPos);
  NTpad = NT(idxPad);
  DVpad = DV(idxPad);

  idx = find(okPos);
  NT = NT(idx);
  DV = DV(idx);
  
  mask = okPos;
  
end