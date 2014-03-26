% Dont use kLow and kHigh for anything!!

function [k,kLow,kHigh] = fitSegregation(obj,distance,segregation,k0)

  if(numel(distance) < 5)
    fprintf('fitSegregation: Warning only has %d data points.\n', numel(distance))
  end

  if(~exist('k0') | numel(k0) ~= 2)
    k0 = [0.5 2];
  end
  
  % opt = statset('RobustWgtFun','bisquare');
  % [k,res,Jac] = nlinfit(distance,segregation,@logistic,k0,opt);

  [k,res,Jac] = nlinfit(distance,segregation,@logistic,k0);
  
  if(0)
    % Debug plots
    x = linspace(min(distance),max(distance),100);
    plot(distance,segregation,'r*', ...
         x,logistic(k,x),'k-');
    pause
  end
  
  CI = nlparci(k,res,'jacobian',Jac);
    
  kLow = CI(:,1);
  kHigh = CI(:,2);
  
  if(all(k == k0))
    disp('k unchanged, k0 was probably VERY poorly choosen')
    keyboard
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function y = logistic(k,x)
    % y = -0.5./(1+(x/inflection).^slope)+1;
    y = 1 - 0.5./(1+(x/k(1)).^k(2));
    
    % If k(1) is negative, we can get imaginary solutions
    % If k(2) is negative, x=0 --> y = 1.
    if(any(k < 0))
      y = ones(size(x))*1e5;
    end
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
end