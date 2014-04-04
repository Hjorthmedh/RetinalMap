% This script sets up the parameter scan for the Grimbert and Cang
% model

runList = 'experiments/GrimbertCang-exploration/runList-%d.txt';
parentCFG = 'experiments/GrimbertCang-wt.txt';
outCFGmask = 'experiments/GrimbertCang-exploration/GrimbertCang-wt-%d.txt';

alphaRange = 0:0.2:1;
betaRange = 1;
thetaRange = 0.1:0.05:0.3;
rhoRange = 0.01:0.02:0.1;

simID = 1;

nJobs = 6;

for i = 1:nJobs
  listName = sprintf(runList, i);
  fidList(i) = fopen(listName,'w');
end

for iA = 1:numel(alphaRange)
  for iB = 1:numel(betaRange)
    for iT = 1:numel(thetaRange)
      for iR = 1:numel(rhoRange)

        fid = fopen(parentCFG,'r');
        fName = sprintf(outCFGmask,simID);
        fidOut = fopen(fName,'w');
        
        str = fgets(fid);
        while(str ~= -1)
          % Copy the string to the output
          fprintf(fidOut, str);
        
          if(findstr(str,'simName'))
            % Add ID to config file name
            str = strrep(str,';',sprintf('-%d;', simID));
          end
          
          str = fgets(fid);
        end
        
        fclose(fid);
        
        % Add extra parameters
        fprintf(fidOut,'obj.info.params.alpha = %f;\n', alphaRange(iA));
        fprintf(fidOut,'obj.info.params.beta = %f;\n', betaRange(iB));        
        fprintf(fidOut,'obj.info.params.theta = %f;\n', thetaRange(iT));
        fprintf(fidOut,'obj.info.params.rho = %f;\n', rhoRange(iR));

        fclose(fidOut);

        jobID = mod(simID-1,nJobs)+1;
        fprintf(fidList(jobID),'runGrimbertCang([],true,1,false,%s);\n',fName); 
        
        simID = simID + 1;
 
        
      end
    end
  end
end

for i = 1:numel(fidList)
  fclose(fidList(i));
end

