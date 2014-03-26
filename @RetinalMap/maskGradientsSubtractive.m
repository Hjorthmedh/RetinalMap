function maskGradientsSubtractive(obj)

  % !!! We are only masking the A system, too little information on
  % B system as it is.
  
  
  % Make sure we are only using the forward system
  assert(obj.typeFlag == 1); 

  % How large a fraction of the reverse gradient should we subtract
  % from the forward gradients?
  %
  % Tsigankov and Koulakov 2006 sets the masking so that it is as
  % large as possible, without there being any overmasking in WT
  
  RGCmaskingFactor = 0.5 % 0.36 gives no overmasking
  SCmaskingFactor = 0.06; % 0.06 gives no overmasking

  % This subtracts the reverse gradients from the forward gradients

  obj.RGCEphA = max(0,obj.RGCEphA - RGCmaskingFactor * obj.RGCephrinA);
  obj.RGCephrinA = max(0,obj.RGCephrinA * (1 - RGCmaskingFactor));  
  
  % obj.RGCEphB = max(0,obj.RGCEphB - RGCmaskingFactor * obj.RGCephrinB);
  % obj.RGCephrinB = max(0,obj.RGCephrinB * (1 - RGCmaskingFactor));


  obj.SCephrinA = max(0,obj.SCephrinA - SCmaskingFactor*obj.SCEphA);
  obj.SCEphA = max(0,obj.SCEphA * (1 - SCmaskingFactor));
 
  % obj.SCephrinB = max(0,obj.SCephrinB - SCmaskingFactor*obj.SCEphB);
  % obj.SCEphB = max(0,obj.SCEphB * (1 - SCmaskingFactor));
  
end