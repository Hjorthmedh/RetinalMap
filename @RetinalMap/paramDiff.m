function allEqual = paramDiff(obj,otherObj)

  % This function ignores all non-char fields of length greater than 1

  allFields = fieldnames(obj); 
  allEqual = true;
  
  for i = 1:numel(allFields)
    fieldA = eval(sprintf('obj.%s',allFields{i}));
    fieldB = eval(sprintf('otherObj.%s',allFields{i}));
    
    if(ischar(fieldA))
      if(~strcmp(fieldA,fieldB))
        % Different strings 
        fprintf('%s : "%s" "%s"\n', allFields{i},fieldA,fieldB)
        allEqual = false;
      end
    elseif(isstruct(fieldA))
      fprintf('Field %s is a struct, ignoring for comparison\n', ...
              allFields{i})
    elseif(numel(fieldA) == 1 & fieldA ~= fieldB)
      fprintf('%s : %f %f\n', allFields{i}, fieldA, fieldB)
      allEqual = false;
    end
  end
      
  if(allEqual)
    disp('paramDiff: all equal')
  end
  
end