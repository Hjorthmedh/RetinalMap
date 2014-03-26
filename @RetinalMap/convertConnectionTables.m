function convertConnectionTables(obj,action)

  switch(action)
    
    case 'post2pre'
      
      % Make sure the post synaptic representation exists
      assert(~isempty(obj.postsynapticConnections))
      
      nConMax = max(max(obj.numPresynapticConnections),obj.maxConnections);
      
      obj.presynapticConnections = zeros(nConMax,obj.nSC,'int32');
      obj.presynapticWeight = zeros(nConMax,obj.nSC);
      
      numCon = zeros(1,obj.nSC);
      
      for i = 1:obj.nRGC
        for j = 1:obj.numPostsynapticConnections(i)
      
          scIdx = obj.postsynapticConnections(j,i);
          scW = obj.postsynapticWeight(j,i);
          
          numCon(scIdx) = numCon(scIdx) + 1;
          
          obj.presynapticConnections(numCon(scIdx),scIdx) = i;
          obj.presynapticWeight(numCon(scIdx),scIdx) = scW;
          
        end
      end
      
      obj.numPresynapticConnections = sum(obj.presynapticWeight > 0,1);
      
    case 'pre2post'
      
      nConMax = max(max(obj.numPostsynapticConnections),obj.maxConnections);
      
      obj.postsynapticConnections = zeros(nConMax,obj.nRGC,'int32');
      obj.postsynapticWeight = zeros(nConMax,obj.nRGC);
    
      numCon = zeros(obj.nRGC,1);
      
      for i = 1:obj.nSC
        for j = 1:obj.numPresynapticConnections(i)
          rIdx = obj.presynapticConnections(j,i);
          rW = obj.presynapticWeight(j,i);
          
          numCon(rIdx) = numCon(rIdx) + 1;
          
          obj.postsynapticConnections(numCon(rIdx),rIdx) = i;
          obj.postsynapticWeight(numCon(rIdx),rIdx) = rW;
          
        end
      end
  
      obj.numPostsynapticConnections = numCon;
            
    case 'mat2pre'
 
      % We need to check how large max connections has to be
      maxCon = max(sum(obj.connectionMatrix > 0,2));

      if(maxCon > obj.maxConnections | isnan(obj.maxConnections))
        fprintf('Setting obj.maxConnections to %d\n', maxCon)
        obj.maxConnections = maxCon;
      end
      
      obj.presynapticConnections = zeros(obj.maxConnections,obj.nSC,'int32');
      obj.presynapticWeight = zeros(obj.maxConnections,obj.nSC);
 
      for i = 1:obj.nSC
        
        % Check that there are no negative values
        assert(nnz(obj.connectionMatrix < 0) == 0);
        
        RGCidx = find(obj.connectionMatrix(i,:) > 0);
        RGCw = obj.connectionMatrix(i,RGCidx);
        
        n = numel(RGCidx);
        obj.presynapticConnections(1:n,i) = RGCidx;
        obj.presynapticWeight(1:n,i) = RGCw;
        
      end
 
      obj.numPresynapticConnections = sum(obj.connectionMatrix > 0,2);
      
      % Need to recalculate totalWeightRGC
      obj.totalWeightRGC = transpose(sum(obj.connectionMatrix,1));
      obj.totalWeightSC = sum(obj.connectionMatrix,2);
      
    case 'pre2mat'
      
      obj.connectionMatrix = sparse(obj.nSC,obj.nRGC);
      
      for i = 1:obj.nSC
        for j = 1:obj.numPresynapticConnections(i)

          RGCidx = obj.presynapticConnections(j,i);
          obj.connectionMatrix(i,RGCidx) = obj.presynapticWeight(j,i);
          
        end        
      end
     
    otherwise
     
      fprintf('Unknown action %s\n (use: pre2post, post2pre, mat2pre', action)
      keyboard
      
  end

end