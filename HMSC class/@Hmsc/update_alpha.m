function alpha = update_alpha(eta,alpha,nf,pi,spatial,spatialMethod,iWgA,iWgRA,detWgA,idDgA,idDW12gA,iFgA,detDgA,alphapwA,an)
% function alpha = update_alpha(eta,alpha,nf,pi,spatial,iWg,iWgR,detWg,Cg,taug,iICCg,alphapw)

% eta - cell array of latent factors at different levels. Each cell corresponds to one level. Each cell contains one
% matrix with dimention of number_of_units_for_this_level to number_of_factors_on this level
% alpha - cell array of parameters of spatial autocorrelation on different levels. Each cell corresponds to one level.
% Within each cell is either an empty object for non-spatial levels, or an array of indices for spatial level. These
% indices are combined with discrite grid for alpha for that level, given by alphapw{level}(:,1)
% nf - array with number of factors at different levels
% pi - matrix of correspondence of observations to units at different levels
% spatial - binary array, indicating whether the level is spatial or not
% iWg - cell array of precalculated inverse matrices
% detWg - cell array of precalculated determinants
% alphapw - cell array of discrite grids for alpha at different levels and the prior probabilities of those
if min(pi)<0
	np = [];
	nr = 0;
else
	[ny nr]=size(pi);
end
% iterate over the latent factors levels, which are defined as spatial
for r=1:nr
	%user specified: layers
	if spatial(r)==true
		alpha1=alpha{r};
		eta1=eta{r};
		detWg = detWgA{r};
		alphapw = alphapwA{r};
		gridN = size(alphapw, 1);
		
		idDg = idDgA{r};
		idDW12g = idDW12gA{r};
		iFg = iFgA{r};
		detDg = detDgA{r};
		
		switch spatialMethod{r}
			case 'GP'
				iWg = iWgA{r};
				iWgR = iWgRA{r};
				tmpMat = permute(sum(mtimesx(iWgR, eta1).^2), [3,2,1]);
			case 'PGP'
				tmpMat2 = mtimesx(eta1',idDW12g);
				tmpMat3 = mtimesx(tmpMat2, iFg);
				tmpMat4 = mtimesx(tmpMat3, permute(tmpMat2,[2,1,3]));
			case 'NNGP'
				iWgR = iWgRA{r};
				tmpMat = nan(gridN,nf(r));
				for ag = 1:gridN
					tmpMat(ag,:) = sum((iWgR{ag}*eta1).^2);
				end
		end
		for h=1:nf(r)
			switch spatialMethod{r}
				case {'GP','NNGP'}
					tmp = tmpMat(:,h);
					like = log(alphapw(:,2))-0.5*detWg-0.5*tmp;
				case 'PGP'
					tmp = nan(gridN, 1);
					
					for ag = 1:gridN
						if alphapw(ag,1) <  1e-5
							tmp(ag) = (eta1(:,h))'*eta1(:,h);
							%detDg1(ag) = 0;
						else
							%idD1 = idD1g1(:,ag);
							%idD1W12 = idD1W12g1(:,:,ag);
							%iF = iFg1(:,:,ag);
							tmp1 = (eta1(:,h))'*(idDg(:,ag).*eta1(:,h));
							%tmp2 = tmpMat2(h,:,ag);
							%tmp2 = (eta1(:,h))'*idD1W12g1(:,:,ag);
							%tmp(ag) = tmp1-tmp2*iFg1(:,:,ag)*tmp2';
							tmp(ag) = tmp1 - tmpMat4(h,h,ag);
						end
					end
					like = log(alphapw(:,2))-0.5*detDg-0.5*tmp;
			end
			like=like-max(like);
			like = exp(like);
			like = like/sum(like);
			alpha1(h)=randsample(length(like),1,true,like);
		end
		alpha{r}=alpha1;
	end
end




%                   C = Cg1(:,:,ag);
%                   tau = taug1(ag);
% %                   M = eye(size(C,2)) + C'*C/tau;
%                   iICC = iICC1(:,:,ag);
%                   d = C'*eta1(:,h);
% %                   tmp(ag) = (eta1(:,h))'*eta1(:,h)/tau - (d'*iICC*d)/tau^2;
%                   tmp(ag) = (eta1(:,h))'*invChol_mex(C*C'+eye(size(C,1))*tau)*eta1(:,h);
