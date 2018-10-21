function etaPred = predictLatentFactor(m, n, r, piCell, xyCell, expected)
calcDirect = false;

res = cell(1, n);
if n==m.postSamN
   postSamInd = 1:n;
else
   postSamInd = randi([1, m.postSamN],1,n);
end

piKey = m.piKey{r};
keySet = unique(piCell);
exInd = m.piMap{r}.isKey(keySet);
if any(~exInd)
   keySetN = keySet(~exInd);
   piKeyN = [piKey; keySetN];
   newMap = containers.Map(keySetN, (length(piKey)+1):length(piKeyN) );
   piMapN = [m.piMap{r}; newMap];
else
   piKeyN = piKey;
   piMapN =  m.piMap{r};
end
piN = cell2mat(piMapN.values(piCell) );
if m.spatial(r)
   piKey = piKeyN;
   xyKey = xyCell(:,1);
   xyMap = containers.Map(xyKey,1:length(xyKey));
   if ~all(xyMap.isKey(piKey))
      error('HMSC: some units, defined at level %d were not given spatial coordinates\n', r);
   end
   ind = cell2mat( xyMap.values(piKey) );
   xy1 = xyCell(ind, 2:size(xyCell, 2));
   xy1 = cell2mat(xy1);
   xy{r} = xy1;
end
pi = piN;

alphaAllPost = cat(1,m.postSamVec(postSamInd).alpha);

newPiInd = ~ismember(pi, m.pi(:,r));
pi1 = pi(newPiInd);
newPiN = length( unique(pi1) );
if m.spatial(r) && newPiN > 0
   if calcDirect == true
      daa = zeros(m.np(r)+newPiN);
      xy1 = xy;
      for j = 1:m.spatDim(r)
         xx = repmat(xy1(:,j), 1, m.np(r)+newPiN);
         dx = xx-xx';
         daa = daa+dx.^2;
      end
      daa = sqrt(daa);
   else
      xys = m.xyKnot{r};
      dss = zeros(size(xys,1));
      for j = 1:m.spatDim(r)
         xx = repmat(xys(:,j), 1, size(xys,1));
         dx = xx-xx';
         dss = dss+dx.^2;
      end
      das = zeros(size(xy1,1),size(xys,1));
      for j = 1:m.spatDim(r)
         xx2 = repmat(xys(:,j), 1, size(xy1,1));
         xx1 = repmat(xy1(:,j), 1, size(xys,1));
         dx = xx1-xx2';
         das = das+dx.^2;
      end
      dss = sqrt(dss);
      das = sqrt(das);
      dns = das(m.np(r)+(1:newPiN),:);
      
      LiFg = cholx_mex(m.iFg{r});
      Wssg = cell(1,size(m.alphapw{r},1));
      iWssg = cell(1,size(m.alphapw{r},1));
      Wnsg = cell(1,size(m.alphapw{r},1));
      WnsiWssg = cell(1,size(m.alphapw{r},1));
      dDng = cell(1,size(m.alphapw{r},1));
      alphaIndUn = sort(unique(cell2mat(alphaAllPost(:,r))));
      for i = 1:length(alphaIndUn)
         ag = alphaIndUn(i);
         alpha = m.alphapw{r}(ag,1);
         if alpha > 0
            Wssg{ag} = exp(-dss/alpha);
            iWssg{ag} = invChol_mex(Wssg{ag});
            Wnsg{ag} = exp(-dns/alpha);
            WnsiWssg{ag} = Wnsg{ag}*iWssg{ag};
            dDng{ag} = ones(size(Wnsg{ag},1),1) - sum(WnsiWssg{ag}.*Wnsg{ag},2);
         else
            Wssg{ag} = eye(size(dss));
            iWssg{ag} = eye(size(dss));
         end
      end
   end
end


for rN = 1:n
   if mod(rN, 100) == 0
      fprintf('Calculating prediction %d\n', rN);
   end
   p = m.postSamVec(postSamInd(rN));
      
   etaM = p.eta{r};
   newPiInd = ~ismember(pi(:,r), m.pi(:,r));
   pi1 = pi(newPiInd, r);
   newPiN = length( unique(pi1) );
   if m.spatial(r) && newPiN > 0
      alphaInd = p.alpha{r};
      alphapw = m.alphapw{r};
      etaN = zeros(newPiN, p.nf(r));
      for j = 1:p.nf(r)
         ag = alphaInd(j);
         alpha = alphapw(ag, 1);
         if alpha > 0
            if calcDirect
               don = daa(1:m.np(r), m.np(r)+(1:newPiN));
               dnn = daa(m.np(r)+(1:newPiN), m.np(r)+(1:newPiN));
               Won = exp(-don/alpha);
               Wnn = exp(-dnn/alpha);
               iWoo = m.iWg{r}(:,:,ag);
               muN = Won' * iWoo * etaM(:,j);
               WN = Wnn - Won' * iWoo * Won;
               WN = (WN+WN') / 2;
               etaN(:,j) = mvnrnd(muN, WN);
               
               LWN = LWNgA{r}{ag};
               etaN(:,j) = muN + LWN*randn(newPiN,1);
            else
               Wss = Wssg{ag};
               Wns = Wnsg{ag};
               WnsiWss = WnsiWssg{ag};
               dDn = dDng{ag};
               iF = m.iFg{r}(:,:,ag);
               idDW12 = m.idDW12g{r}(:,:,ag);
               LiF = LiFg(:,:,ag);
               
               muS1 = iF*idDW12'*etaM(:,j);
               epsS1 = LiF*randn(size(Wss,1),1);
%                etaS = Wss*(muS1+(~expected)*epsS1);
%                muN = WnsiWss*etaS;
               muN = Wns*(muS1+(~expected)*epsS1);
               etaN(:,j) = muN + (~expected)*sqrt(dDn).*randn(newPiN,1);               
            end
         else
            etaN(:,j) = (~expected)*normrnd(0,1,[newPiN,1]);
         end
      end
   else
      etaN = (~expected)*normrnd(0,1,[newPiN,p.nf(r)]);
   end
   eta = [etaM; etaN];
   res{rN} = eta(pi,:);
end
etaPred = res;

end


