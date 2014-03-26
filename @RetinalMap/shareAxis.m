% This function makes all axis for multiple figures the same

function shareAxis(figs)
  
  minX = inf;
  maxX = -inf;
  minY = inf;
  maxY = -inf;
  
  for i = 1:numel(figs)
    a = axis(figs(i));
    minX = min(minX,a(1));
    maxX = max(maxX,a(2));
    minY = min(minY,a(3));
    maxY = max(maxY,a(4));
  end

  for i = 1:numel(figs)
    figure(figs(i))
    axis([minX maxX minY maxY]);
  end
end