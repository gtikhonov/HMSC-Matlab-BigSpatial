function [ Ypred ] = crossValidate(m,part,mcmcRepN,postRep,postThin,predN)
partUn = sort(unique(part));
partN = length(partUn);
Ypred = nan(size(m.Y));

mC = m.copy();
mC.postRemove();
for pN = 1:partN
   fprintf('Crossvalidation - partition %d out of %d\n', pN, partN);
   indVal = (part==partUn(pN));
   indTrain = ~indVal;
%    mT = Hmsc(m.folder, m.traits, m.speciesX, m.phylogeny, m.includeXs, m.includeXv, m.outlierSpecies, m.spatial, m.factorCov);
   mT = mC.copy();
   YT = m.Y(indTrain,:);
   piCellT = m.piCell(indTrain,:);
   YV = m.Y(indTrain,:);
   piCellV = m.piCell(indVal,:);
   if m.speciesX == false
      XT = m.X(indTrain,:);
      XV = m.X(indVal,:);
   else
      XT = cell(1,m.ns);
      XV = cell(1,m.ns);
      for j = 1:m.ns
         XT{j} = m.X{j}(indTrain,:);
         XV{j} = m.X{j}(indVal,:);
      end
   end
   if m.includeXs
      XsT = m.Xs(indTrain,:);
      XsV = m.Xs(indVal,:);
   else
      XsT = [];
      XsV = [];
   end
   if m.includeXv
      XvT = m.Xv(indTrain,:);
      XvV = m.Xv(indVal,:);
   else
      XvT = [];
      XvV = [];
   end
   if m.traits
      TT = m.T;
   else
      TT = [];
   end
   
   mT.setData(YT,m.dist(:,1:2),XT,piCellT,m.xyCell,TT,m.C,m.XrCell,XsT,XvT);
   mT.sampleMCMC(mcmcRepN, false, [], 2);
   mT.setPostThinning(postRep, postThin);
   
   predTestCell = m.predict(predN, XV, XsV, XvV, piCellV, m.xyCell, m.XrCell, true);
   Ypred(indVal,:) = mean(cat(3, predTestCell{:}),3);
end

end

