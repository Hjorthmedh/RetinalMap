function d = SCDistance(obj)

  memoryFlag = 1; % Tries to conserve memory
  
  if(memoryFlag)
    
    d = zeros(numel(obj.SCap),numel(obj.SCap));
    
    for i = 1:numel(obj.SCap)
      d(:,i) = sqrt((obj.SCap-obj.SCap(i)).^2 + (obj.SCml-obj.SCml(i)).^2);
    end    
    
  else
  
    d = sqrt((kron(obj.SCap,ones(1,length(obj.SCap))) ...
              - kron(ones(length(obj.SCap),1),transpose(obj.SCap))).^2 ...
             + (kron(obj.SCml,ones(1,length(obj.SCml))) ...
                - kron(ones(length(obj.SCml),1), ...
                       transpose(obj.SCml))).^2);
    
  end
  
end
