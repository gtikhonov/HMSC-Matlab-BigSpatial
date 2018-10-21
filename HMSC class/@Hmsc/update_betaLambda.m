function [betaNew,lambdaNew] = update_betaLambda(Y,X,Xr,Xs,Xv,z,beta,nf,pi,ncr,eta,lambda,psijh,delta,etas,lambdas,betav,qv,gamma,sigma,T,ph,iV,rho,phylogeny,iQg,detQg,outlierspecies,spatial,factorcov,speciesX,includeXs,includeXv,an)

[ny, ns] = size(Y);
isigmaDiag = diag(sigma).^-1;
Ys = (Y==1) - (Y==0);
if speciesX
   nc=size(X{1},2);
else
   nc=size(X,2);
end
if isempty(spatial)
   nr = 0;
else
   [~, nr]=size(pi);
end
Sr = zeros(ny,ns);
if includeXs == 1
   Sr = Sr + Xs*etas*lambdas;
elseif includeXs == 2
   for j=1:ns
      Sr(:,j) = Sr(:,j) + Xs{j}*etas*lambdas(:,j);
   end
end
if includeXv == 1
   Sr = Sr + Xv*(qv.*betav);
end
S = z - Sr;

nfCumSum = [0,cumsum(nf)];
nfSum = nfCumSum(end);
if ~phylogeny
   XEta = nan(ny,nc+nfSum);
   XEta(:,1:nc) = X;
   for r = 1:nr
      XEta(:,(nc+nfCumSum(r)+1):(nc+nfCumSum(r+1))) = eta{r}(pi(:,r),:);
   end
   zRand = randn([nc+nfSum,ns]);
   gaT = gamma'*T';
   priorMean = [gaT;zeros(nfSum,ns)];
   priorPrecMat = zeros(nc+nfSum,nc+nfSum,ns);
   priorPrecMat(1:nc,1:nc,:) = repmat(iV,[1,1,ns]);
   for r = 1:nr
      priorLambda = bsxfun(@times,psijh{r},cumprod(delta{r})');
      ind1 = repmat((nc+nfCumSum(r)+(1:nf(r))-1)*(nc+nfSum)+(nc+nfCumSum(r)+(1:nf(r))),[1,ns]);
      ind2 = repmat(((nc+nfSum)^2*((1:ns)-1)),[nf(r),1]);
      ind = ind1 + ind2(:)';
      priorLambda = priorLambda';
      priorPrecMat(ind) = priorLambda(:);
   end
   
   XEta2 = XEta'*XEta;
   iSigmaXEta2 = repmat(permute(isigmaDiag,[2,3,1]),[nc+nfSum,nc+nfSum]).*repmat(XEta2,[1,1,ns]);
   RiU = permute(cholx_mex(priorPrecMat + iSigmaXEta2),[2,1,3]);
   SbR = S;
   
   RU = invTrix_mex(RiU,'U');
   U = mtimesx(RU,permute(RU,[2,1,3]));
   m = mtimesx(U, mtimesx(priorPrecMat,permute(priorMean,[1,3,2])) + permute(repmat(isigmaDiag',nc+nfSum,1).*(XEta'*SbR),[1,3,2]));
   betaLambdaNew = permute(m + mtimesx(RU,permute(zRand,[1,3,2]))*sqrt(an), [1,3,2]);
   % 	fprintf('%e\n', max(max(abs(betaLambdaNew - betaLambdaNew1))));
   
   betaNew = betaLambdaNew(1:nc,:);
   lambdaNew = mat2cell(betaLambdaNew(nc+1:end,:), nf)';
else
   error('HMSC: update_betaLambda is not yet written for phylogenic models');
end
end


% species by species samplimng
%       SbR = (S - bR{1})./bR{2};
%    betaLambdaNew = nan(nc+nfSum,ns);
%    for j = 1:ns
%       priorMean = (T(j,:)*gamma)';
%       priorMean = [priorMean; zeros(nfSum,1)];
%       priorPrecMat = zeros(nc+nfSum);
%       priorPrecMat(1:nc,1:nc) = iV;
%       for r = 1:nr
%          priorLambda = bsxfun(@times,psijh{r},cumprod(delta{r})');
%          priorPrecMat(nc+nfCumSum(r)+1:nc+nfCumSum(r+1), nc+nfCumSum(r)+1:nc+nfCumSum(r+1)) = diag(priorLambda(j,:));
%       end
%       sqrtRXEta = XEta ./ repmat(sqrt(bR{2}(:,j)),[1,nc+nfSum]);
%       RXEta2 = sqrtRXEta'*sqrtRXEta;
%       iSigmaRXEta2 = isigmaDiag(j)*RXEta2;
%       iU = priorPrecMat + iSigmaRXEta2;
%       LiU = chol(iU,'lower');
%       tmp1 = priorPrecMat*priorMean + XEta'*SbR(:,j);
%       tmp2 = LiU\tmp1;
%       betaLambdaNew(:,j) = LiU'\(tmp2+zRand(:,j));
%    end
%    betaLambdaNew1 = betaLambdaNew;