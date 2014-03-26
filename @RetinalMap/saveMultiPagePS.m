% For first page, set appendFlag = false, for subsequent pages set
% it to true

function saveMultiPagePS(obj,fileName,fig,appendFlag,format)

  if(~exist('appendFlag'))
    appendFlag = 0;
  end
  
  if(~exist('format'))
    format = 'eps';
  end

  switch(format)
    case 'eps'
      formatStr = '-dpsc2';
    case 'pdf'
      formatStr = '-dpdf';
    otherwise
      fprintf('Unknown format %s, defaulting to pdf\n', format)
      formatStr = '-dpdf';
  end
    
  
  for i = 1:numel(fig)
    if(i == 1 & ~appendFlag)
      print(fig(i),formatStr, fileName)      
    else
      print(fig(i),formatStr, '-append', fileName)
    end
  end

end