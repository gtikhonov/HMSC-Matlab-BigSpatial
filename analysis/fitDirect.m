restoredefaultpath;
readData;
execStartTime = tic();
RCVec = 1:2;
modelType = 'Direct';

mcmcRepN = 10;
mcmcSam = 100;
mcmcThin = 100;
nf = 2;

addpath(fullfile('..', 'HMSC class'));
baseFolder = '.';

for RC = RCVec
    for spN = 1:size(speciesPartitionUsed,2)
        for sN = 1:size(dataPartitionUsed,2)
            if speciesPartitionUsed(2,spN)==0 || dataPartitionUsed(2,sN)==0 || dataPartitionUsed(1,sN)>1600
                continue;
            end
            rng(RC)
            fitSize = dataPartitionUsed(1,sN);
            fitInd = dataPartition>0 & dataPartition<=sN;
            X = [ones(fitSize,1), XAll(fitInd,:), XAll(fitInd,:).^2];
            
            spSize = speciesPartitionUsed(1,spN);
            speciesInd = speciesPartition>0 & speciesPartition<=spN;
            Y = YAll(fitInd,speciesInd);
            piCell = piCellAll(fitInd,:);
            
            addPath1 = sprintf('%s %.2d',modelType,nf);
            addPath2 = sprintf('%.4d-%.5d',spSize,fitSize);
            addPath3 = sprintf('RC %d, repN %.5d, sam %.4d, thin %.3d', RC, mcmcRepN, mcmcSam, mcmcThin);
            addPath = fullfile(addPath1, addPath2, addPath3);
            folder = fullfile(baseFolder, 'fitted models',  addPath);
            fprintf('%s\n', addPath);
            if exist(fullfile(folder,'fm.mat'),'file')
                continue;
            end
            tempFileName = sprintf('%s --- %s.txt', addPath2, addPath3);
            tempFileFullPath = fullfile(baseFolder,'fitted models',addPath1,tempFileName);
            if exist(tempFileFullPath,'file')
                continue;
            end
            [~, ~, ~] = mkdir(fullfile(baseFolder,'fitted models',addPath1));
            fid = fopen(tempFileFullPath, 'wt');
            fprintf(fid,'%s --- %s', getenv('computername'), datestr(datetime('now','Format','y-MM-dd HH:mm:ss')));
            fclose(fid);
            
            m = Hmsc(folder, true, false, false, 0, false, false, [true], [false]);
            m.setData(Y,'probit',X,piCell,{xyCellAll},T(speciesInd,:),[],[],[],[]);
            m.setCovNames(['Intercept',covNames,strcat(covNames,'_2')]);
            m.setSpeciesNames(spNames(speciesInd));
            m.setLevelNames({'site'});
            m.setTraitNames(traitNames);
            m.setCovScaling([2,ones(1,size(XAll,2)*2)]);
            m.setSpatialMethod({'GP'});
            
            m.setPriorsDefault();
            m.setMCMCOptions(mcmcSam, mcmcThin);
            m.setMCMCAdapt([0,0], nf);
            m.setMCMCSaveOptions(true, false);
            anVec = ones(1,mcmcRepN);
            startTimer = tic;
            m.sampleMCMC(mcmcRepN, false, [], 3, true, anVec);
            execT = toc(startTimer);
            
            [s, mess, messid] = mkdir(folder);
            save(fullfile(folder,'fm.mat'),'m','execT');
            save(fullfile(folder,'execTime.mat'),'execT');
            delete(tempFileFullPath);
        end
    end
end
toc(execStartTime)