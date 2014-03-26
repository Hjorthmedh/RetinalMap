function [AP,ML,APpad,MLpad,mask] = filterSCPosition(obj,AP,ML)

  if(isempty(obj.SCmask))
    obj.SCmask = imread('images/SCadultmouse.png');
    obj.SCmask = obj.SCmask(:,:,2) < 50;
  end
    
  w = size(obj.SCmask,2);
  h = size(obj.SCmask,1);
  
  mask = zeros(size(AP));
  
  imgMask = 1 <= ceil(ML*h) & ceil(ML*h) <= size(obj.SCmask,1) ...
            & 1 <= ceil((1-AP)*w) & ceil((1-AP)*w) <= size(obj.SCmask,2);

  idxIn = find(imgMask);
  
  idx = sub2ind(size(obj.SCmask),ceil(ML(idxIn)*h),ceil((1-AP(idxIn))*w));

  mask(idxIn) = obj.SCmask(idx);  
  
  okIdx = find(mask);
  outsideIdx = find(~mask);
  
  APpad = AP(outsideIdx);
  MLpad = ML(outsideIdx);
  
  AP = AP(okIdx);
  ML = ML(okIdx);
  
end