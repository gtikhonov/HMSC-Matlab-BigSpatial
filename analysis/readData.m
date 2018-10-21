clearvars;
restoredefaultpath;
RC = 1;
rng(RC);
dataFolder = './data';

%%%%% reading data
TX = readtable(fullfile(dataFolder,'X.csv'));
N = size(TX,1);
XAll = table2array(TX);
covNames = TX.Properties.VariableNames;

piCellAll = cellfun(@num2str, num2cell((1:N)'), repmat({'%.5d'},N,1),'UniformOutput', false);

xy = csvread(fullfile(dataFolder,'xy.csv'));
xyCellAll = [piCellAll, num2cell(xy)];

TY = readtable(fullfile(dataFolder,'Y.csv'));
spNames = TY.Properties.VariableNames;
YAll = table2array(TY);
% ind5 = sum(YAll)>=5;

TT = readtable(fullfile(dataFolder,'traits.csv'));
traitNames = TT.Properties.VariableNames;
T = [ones(size(TT,1),1), table2array(TT(:,2:end))];

dataPartition = csvread(fullfile(dataFolder,'data partition.csv'));
dataPartitionUsed = csvread(fullfile(dataFolder,'data partition used.csv'));
speciesPartition = csvread(fullfile(dataFolder,'species partition.csv'));
speciesPartitionUsed = csvread(fullfile(dataFolder,'species partition used.csv'));

% dataPartitionUsed(2,end) = 0;
% speciesPartitionUsed(2,end) = 0;
speciesPartitionUsed(2,1) = 0;

dataIndexLength = sum(dataPartitionUsed(2,:));
speciesIndexLength = sum(speciesPartitionUsed(2,:));

knotsNumber = [16,64,256,1024];
nnNumber = [10,20];
hull = csvread(fullfile(dataFolder,'hull.csv'));

 
% spNames = cellfun(@num2str, num2cell((1:size(TY,2))'), repmat({'sp%.4d'},size(TY,2),1),'UniformOutput', false);
% TY.Properties.VariableNames = spNames;
% TT.CommonName = spNames;
% writetable(TY,fullfile(dataFolder,'Y1.csv'));
% writetable(TT,fullfile(dataFolder,'T1.csv'));
% readtable(fullfile(dataFolder,'Y.csv'));