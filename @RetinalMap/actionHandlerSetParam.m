function actionHandlerSetParam(obj,aParam)

  fprintf('actionHandlerSetParam called at epoch %d\n', obj.curStep/obj.nSC)

  try
    for i = 1:numel(aParam)
      varName = aParam{i}{1};
      varValue = aParam{i}{2};
      cmd = sprintf('obj.%s = varValue;',varName);
      eval(cmd);
    end
  catch e
    getReport(e)
    keyboard
  end
  
end