function eta = update_eta(z,X,Xr,Xs,Xv,beta,sigma,eta,alpha,lambda,etas,lambdas,betav,qv,nf,pi,ncr,spatial,spatialMethod,factorcov,iWgA,idDgA,idD1W12gA,FgA,speciesX,includeXs,includeXv,an)

[ny ns] = size(z);
isigma = diag(diag(sigma).^-1);
if length(spatial)==0
	np = [];
	nr = 0;
else
	[ny nr]=size(pi);
	np = max(pi);
end
for r=1:nr
	S = z;
	if speciesX
		for j=1:ns
			S(:,j) = S(:,j) - X{j}*beta(:,j);
		end
	else
		S = S - X*beta;
	end
	if includeXs == 1
		S = S - Xs*etas*lambdas;
	elseif includeXs == 2
		for j=1:ns
			S(:,j) = S(:,j) - Xs{j}*etas*lambdas(:,j);
		end
	end
	if includeXv == 1
		S = S - Xv*(qv.*betav);
	end
	
	for r2=1:nr
		if ~(r==r2)
			eta1=eta{r2};
			lambda1=lambda{r2};
			if factorcov(r2) > 0
				Xr1 = Xr{r2};
				for k=1:ncr(r2)
					if factorcov(r2) == 1
						XrEta = repmat(Xr1(:,k), 1, nf(r2)).*eta1;
						S = S-XrEta(pi(:,r2),:)*lambda1(:,:,k);
					else
						for j=1:ns
							XrEta = repmat(Xr1(:,k), 1, nf(r2)).*eta1;
							S(:,j) = S(:,j) - XrEta(pi(:,r2),:)*lambda1(:,j,k);
						end
					end
				end
			else
				S = S-eta1(pi(:,r2),:)*lambda1;
			end
		end
	end
	lpi = pi(:,r);
	unLpi = unique(lpi);
	lambda1=lambda{r};
	if factorcov(r)
		Xr1 = Xr{r};
		if ~spatial(r)
			eta1=eta{r};
			if length(unLpi)==ny % i1 corresponds to residuals, use a faster updater
				for q = 1:ny
					p = lpi(q);
					lambdaX = zeros(nf(r), ns);
					for k = 1:ncr(r)
						if factorcov(r) == 1
							lambdaX = lambdaX + Xr1(p,k)*lambda1(:,:,k);
						else
							for j=1:ns
								lambdaX(:,j) = lambdaX(:,j) + Xr1{j}(p,k)*lambda1(:,j,k);
							end
						end
					end
					lambda2 = (lambdaX*isigma*lambdaX');
					Veta1 = eye(nf(r)) + lambda2;
					Tx = cholcov(Veta1); [~,R] = qr(Tx);
					iR = inv(R); Veta = iR*iR';  % Veta = inv(Veta1)
					Meta = S(q,:)*isigma*lambdaX'*Veta;
					eta1(p,:) = Meta + randn([1,nf(r)])*iR';
				end
			else
				for q = 1:length(unLpi)
					p = unLpi(q);
					lambdaX = zeros(nf(r), ns);
					for k = 1:ncr(r)
						if factorcov(r) == 1
							lambdaX = lambdaX + Xr1(p,k)*lambda1(:,:,k);
						else
							for j=1:ns
								lambdaX(:,j) = lambdaX(:,j) + Xr1{j}(p,k)*lambda1(:,j,k);
							end
						end
					end
					lambdaXT = lambdaX';
					lambda2 = (lambdaX*isigma*lambdaXT);
					rows = lpi==p;
					Veta1 = eye(nf(r)) + lambda2*sum(rows);
					Tx = cholcov(Veta1);
					[~,R] = qr(Tx);
					iR = inv(R);
					iRT = iR';
					Veta = iR*iRT;  % Veta = inv(Veta1)
					Meta = sum(S(rows,:),1)*isigma*lambdaXT*Veta;
					%Meta = sum(S(rows,:),1)*isigma*lambdaX'/Veta1;
					eta1(p,:) = Meta + randn([1,nf(r)])*iRT;
				end
			end
		else % SPATIAL LATENT FACTORS
			iWg = iWgA{r};
			alpha1=alpha{r};
			iW = zeros( np(r)*nf(r) );
			for h=1:nf(r)
				iW((h-1)*np(r)+(1:np(r)), (h-1)*np(r)+(1:np(r))) = iWg(:,:,alpha1(h));
			end
			Xrk = cell(1, ncr(r));
			isigmaLambdak = cell(1, ncr(r));
			for k = 1:ncr(r)
				if factorcov(r) == 1
					Xrk1 = diag(Xr1(:,k));
					Xrk{k} = Xrk1(lpi,:);
				else
					fprintf('Not implemented yet\n');
					Xrk1 = cell(1, ns);
					for j = 1:ns
						Xrk1{j} = diag(Xr1{j}(:,k));
						Xrk1{j} = Xrk1{j}(lpi,:);
					end
					Xrk{k} = Xrk1;
				end
				isigmaLambdak{k} = isigma*lambda1(:,:,k)';
			end
			Ueta1 = iW;
			for k1 = 1:ncr(r)
				Xrk1 = Xrk{k1};
				for k2 = 1:ncr(r)
					Xrk2 = Xrk{k2};
					if factorcov(r) == 1
						Ueta1 = Ueta1 + kron(lambda1(:,:,k1)*isigmaLambdak{k2}, Xrk1'*Xrk2);
					else
						fprintf('Not implemented yet\n');
					end
				end
			end
			Ueta = inv(Ueta1);
			Ueta = (Ueta+Ueta')/2;
			Meta = zeros(nf(r)*np(r), 1);
			for k1 = 1:ncr(r)
				Xrk1 = Xrk{k1};
				vec = Xrk1' * S * isigmaLambdak{k1};
				Meta = Meta + Ueta*vec(:);
			end;
			feta = mvnrnd(Meta,Ueta);
			eta1 = reshape(feta,[np(r),nf(r)]);
		end
	else
		if spatial(r)==0 % NON SPATIAL LATENT FACTORS
			eta1=eta{r};
			if length(unLpi)==ny % i1 corresponds to residuals, use a faster updater
				Veta1 = eye(nf(r)) + lambda1*isigma*lambda1';
				Tx = cholcov(Veta1);
				[~,R] = qr(Tx);
				iR = inv(R);
				Veta = iR*iR';
				Meta = S*isigma;
				Meta = Meta*lambda1';
				Meta = Meta*Veta;                        % ny x k
				tmp = Meta + randn([ny,nf(r)])*iR'*sqrt(an);       % update eta1 in a block
				eta1(lpi,:) = tmp;
			else % i1 corresponds to a random effect, use a more general updater
				lambda2 = (lambda1*isigma*lambda1');
				for q = 1:length(unLpi)
					p = unLpi(q);
					rows = lpi==p;
					Veta1 = eye(nf(r)) + lambda2*sum(rows);
					Tx = cholcov(Veta1); 	[~,R] = qr(Tx);
					iR = inv(R); Veta = iR*iR';   % Veta = inv(Veta1)
					Meta = sum(S(rows,:),1)*isigma*lambda1'*Veta;
					eta1(p,:) = Meta + randn([1,nf(r)])*iR';
				end
			end
		else % SPATIAL LATENT FACTORS
			alpha1=alpha{r};
			switch spatialMethod{r}
				case 'GP'
					iWg = iWgA{r};
					iW = zeros( np(r)*nf(r) );
					for h=1:nf(r)
						iW((h-1)*np(r)+(1:np(r)), (h-1)*np(r)+(1:np(r))) = iWg(:,:,alpha1(h));
					end
					% may be try to change to block diagonal construction
					%c=num2cell(A,[1 2])
					%M=blkdiag(c{:})
					if sum(lpi==(1:ny)')==ny % i1 corresponds to spatial residuals, use a faster updater
						iSigmaLambdaT = isigma*lambda1';
						tmp1 = lambda1*iSigmaLambdaT;
						tmp1s = kron(tmp1,speye(ny));
						iUeta = iW + tmp1s;
						fS = S*iSigmaLambdaT;
						L = chol(iUeta, 'lower');
						opts1.LT = true;
						opts2.LT = true; opts2.TRANSA = true;
						tmp2 = linsolve(L,fS(:),opts1) + randn([np(r)*nf(r),1]);
						feta = linsolve(L, tmp2, opts2);
						eta1 = reshape(feta,[np(r),nf(r)]);
					else % i1 corresponds to a spatial random effect, use a more general updater
						Pi = eye(np(r));
						Pi = Pi(lpi, :);
						iSigmaLambdaT = repmat(iSigmaDiag,1,nf(r)).*lambda1';
						tmp1 = lambda1*iSigmaLambdaT;
						tmp1s = kron(tmp1,diag(sum(Pi)));
						iUeta = iWs + tmp1s;
						fS = (Pi'*S)*iSigmaLambdaT;
						L = chol(iUeta, 'lower');
						opts1.LT = true;
						opts2.LT = true; opts2.TRANSA = true;
						tmp2 = linsolve(L,fS(:),opts1) + randn([np(r)*nf(r),1]);
						feta = linsolve(L, tmp2, opts2);
						eta1 = reshape(feta,[np(r),nf(r)]);
					end
				case 'PGP'
					% currently only for spatial residuals
					idDg = idDgA{r};
					idD1W12g = idD1W12gA{r};
					Fg = FgA{r};
					nK = size(Fg,1);
					idD = idDg(:,alpha1);
					
					F = zeros( nK*nf(r) );
					idD1W12 = zeros(np(r)*nf(r), nK*nf(r));
					for h=1:nf(r)
						F((h-1)*nK+(1:nK), (h-1)*nK+(1:nK)) = Fg(:,:,alpha1(h));
						idD1W12((h-1)*np(r)+(1:np(r)), (h-1)*nK+(1:nK)) = idD1W12g(:,:,alpha1(h));
					end
					tmp = isigma*lambda1';
					fS = S*tmp;
					fS = fS(:);
					LamSigLamT = lambda1*tmp;
					
% 					A0 = LamSigLamT + diag(idD);
% 					iA0 = invChol_mex(A0);
% 					iA = kron(iA0, speye(ny));
% 					LiA0 = chol(iA0, 'lower');
% 					LiA = kron(LiA0, speye(ny));
% 					iA1 = inv(kron(LamSigLamT, speye(ny))+diag(idD(:))); % just
% 					check
					
					B0 = repmat(LamSigLamT,[1,1,ny]);
					idDV = reshape(permute(idD,[2,1]), [1,ny*nf]);
					tmp = repmat(((1:ny)-1), nf,1);
					ind = (tmp(:)*nf^2)' + repmat(nf*((1:nf)-1)+(1:nf), [1,ny]);
					B0(ind) = B0(ind) + idDV;
					
					B1 = invCholx_mex(B0);
					tmp1 = repmat(((1:nf)-1),ny,1)*ny;
					ind1 = repmat(repmat(1:ny,1,nf) + tmp1(:)', 1, nf);
					tmp1 = repmat(((1:nf)-1),ny*nf,1)*ny;
					ind2 = repmat(1:ny,1,nf^2) + tmp1(:)';
					iA = sparse(ind1, ind2, reshape(permute(B1, [3,1,2]), ny*nf^2, 1));
					LB1 = cholx_mex(B1);
					LiA = sparse(ind1, ind2, reshape(permute(LB1, [3,1,2]), ny*nf^2, 1));
					
					iAidD1W12 = iA*idD1W12;
					H = F - idD1W12'*iAidD1W12;
					RH = chol(H);
					iRH = invTri_mex(RH,'U');
					% iRH = chol(iH,'lower');
					
					mu1 = iA*fS;
					tmp1 = iAidD1W12*iRH;
					mu2 = tmp1*(tmp1'*fS);
					
					etaR = LiA*randn(np(r)*nf(r),1) + tmp1*randn(nK*nf(r),1);
					eta1 = reshape(mu1+mu2+etaR,[np(r),nf(r)]);
				case 'NNGP'
					iWg = iWgA{r};
					iW = sparse(np(r)*nf(r), np(r)*nf(r));
% 					for h=1:nf(r)
% 						iW((h-1)*np(r)+(1:np(r)), (h-1)*np(r)+(1:np(r))) = iWg{alpha1(h)};
% 					end
% 					tmp = isigma*lambda1';
% 					tmp1 = lambda1*tmp;
% 					tmp1s = kron(tmp1,speye(ny));
% 					iUeta = iW + tmp1s;
					
					for h=1:nf(r)
						iW(h+(0:np(r)-1)*nf(r), h+(0:np(r)-1)*nf(r)) = iWg{alpha1(h)};
					end
					tmp = isigma*lambda1';
					tmp1 = lambda1*tmp;
					tmp1s = kron(speye(ny),tmp1);
					iUeta = (iW + tmp1s);
					
					ran = randn([np(r)*nf(r),1]);
					
% 					LiU = chol(iUeta,'lower');
% 					fS = (S*tmp)';
% 					feta = LiU'\(LiU\fS(:) + ran);
% 					eta1a = reshape(feta,[nf(r),np(r)])';
					
					[LiU,p,q] = chol(iUeta,'lower');
					if p > 0
						warning('Covariance matrix for eta is not positive definite')
					end
					fS = (S*tmp)';
					feta = q*(LiU'\(LiU\(q'*fS(:)) + ran));
					eta1 = reshape(feta,[nf(r),np(r)])';
					
					
% 					nnz(LiU) / numel(LiU)
% 					subplot(1,2,1)
% 					spy(LiU)
% 					[LiU,p,P] = chol(iUeta,'lower');
% 					nnz(LiU) / numel(LiU)
% 					subplot(1,2,2)
% 					spy(LiU)
			end
		end
	end
	eta{r} = eta1;
end
