function maskGradientsMultiplicative(obj)

  % !!! We are only masking the A system, too little information on
  % B system as it is.
  
  
  % Make sure we are only using the forward system
  % assert(obj.typeFlag == 1); 

  % How large a fraction of the reverse gradient should we subtract
  % from the forward gradients?
  %
  % Tsigankov and Koulakov 2006 sets the masking so that it is as
  % large as possible, without there being any overmasking in WT
  
  RGCmaskingFactor = 0.5 
  SCmaskingFactor = 0.5; 

  % This subtracts the reverse gradients from the forward gradients

  % The min&max are to prevent overmasking
  maskAR = min(RGCmaskingFactor*obj.RGCEphA .* obj.RGCephrinA, ...
              max(obj.RGCEphA,obj.RGCephrinA));
  
  obj.RGCEphA = obj.RGCEphA - maskAR;
  obj.RGCephrinA = obj.RGCephrinA - maskAR;
  
  maskAS = min(SCmaskingFactor*obj.SCephrinA.*obj.SCEphA, ...
              max(obj.SCEphA,obj.SCephrinA));

  
  obj.SCephrinA = obj.SCephrinA - maskAS;
  obj.SCEphA = obj.SCEphA - maskAS;
  
end