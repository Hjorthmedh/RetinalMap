function d = retinalDistance(obj)

  memoryFlag = 1; % Minimize the memory footprint
  
  switch(obj.eyeType)
  
    case 'sphere'
      
      x = sin(obj.RGCtheta).*cos(obj.RGCphi);
      y = sin(obj.RGCtheta).*sin(obj.RGCphi);
      z = cos(obj.RGCtheta);

      if(memoryFlag)
      
        d = zeros(numel(x),numel(x));
        
        for i = 1:numel(x)
          d(:,i) = acos(x*x(i) + y*y(i) + z*z(i));
        end
        
        % !!!
        
        d2 = acos(kron(x,ones(1,length(x))) ...
                  .* kron(ones(length(x),1),transpose(x)) ...
                  + kron(y,ones(1,length(y))) ...
                  .* kron(ones(length(y),1),transpose(y)) ...
                  + kron(z,ones(1,length(z))) ...
                  .* kron(ones(length(z),1),transpose(z)));
        
        
        if(nnz(d-d2))
          disp('retinalDistance for sphere incorrect!')
          keyboard
        else
          disp('retinalDistance for sphere correct - remove debug code')
        end
        
      else
      
        % arccos of the scalar product of the two vectors
        d = acos(kron(x,ones(1,length(x))) ...
                 .* kron(ones(length(x),1),transpose(x)) ...
                 + kron(y,ones(1,length(y))) ...
                 .* kron(ones(length(y),1),transpose(y)) ...
                 + kron(z,ones(1,length(z))) ...
                 .* kron(ones(length(z),1),transpose(z)));
      end
        
      % For diagonal elements we sometimes get small imaginary numbers
      % Per definition diagonal should be zero.

      for i = 1:size(d,1)
        d(i,i) = 0;
      end
      
    case 'disk'
  
      if(memoryFlag)
      
        d = zeros(numel(obj.RGCnt),numel(obj.RGCnt));
        
        for i = 1:numel(obj.RGCnt)
          d(:,i) = sqrt((obj.RGCnt-obj.RGCnt(i)).^2 ...
                        + (obj.RGCdv-obj.RGCdv(i)).^2);
        end
        
      else

        d = sqrt((kron(obj.RGCnt,ones(1,length(obj.RGCnt))) ...
                  - kron(ones(length(obj.RGCnt),1), ...
                         transpose(obj.RGCnt))).^2 ...
                 + (kron(obj.RGCdv,ones(1,length(obj.RGCdv))) ...
                    - kron(ones(length(obj.RGCdv),1), ...
                           transpose(obj.RGCdv))).^2);
        
        
      end
      
      
      
     
    otherwise
      fprintf('retinalDistance: Unknown eye type %s\n', obj.eyeType)
      keyboard
  
        end
  
  if(nnz(imag(d)))
    disp('Imaginary distances --- how?')
    keyboard
  end

end
