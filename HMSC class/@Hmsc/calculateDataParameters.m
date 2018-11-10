function [iQg,detQg,iWgA,iWgRA,detWgA,idDgA,idDW12gA,FgA,iFgA,detDgA] = calculateDataParameters(m)
% function [iQg,detQg,iWgA,iWgAR,detWgA,CgA,taugA,iICCgA] = calculateDataParameters(m)

iQg = nan(m.ns,m.ns,size(m.rhopw,1));
detQg = nan(size(m.rhopw,1),1);
iWgA = cell(1,m.nr);
iWgRA = cell(1,m.nr); %cholesky component of iW
detWgA = cell(1,m.nr);
idDgA = cell(1,m.nr);
idDW12gA = cell(1,m.nr);
FgA = cell(1,m.nr);
iFgA = cell(1,m.nr);
detDgA = cell(1,m.nr);

if m.phylogeny
	rhopw = m.rhopw;
	C = m.C;
	iC = corrcov(inv(m.C));
	for rg=1:size(rhopw, 1)
		rho = rhopw(rg, 1);
		if rho >= 0
			rhoC = rho*C;
		else
			rhoC = (-rho)*iC;
		end
		Q = rhoC+(1-abs(rho))*eye(m.ns);
		iQg(:,:,rg) = inv(Q);
		cQ = chol(Q);
		detQg(rg) = 2*sum(log(diag(cQ)));
	end
end

for r=1:m.nr
	if m.spatial(r)
		alphapw=m.alphapw{r};
		np = m.np(r);
		xy = m.xy{r};
		alphagN = size(alphapw, 1);
		
		iWg = [];
		iWgR = [];
		detWg = [];
		idDg = [];
		idDW12g = [];
		Fg = [];
		iFg = [];
		detDg = [];
		
		switch m.spatialMethod{r}
			case 'GP'
				di = zeros(np);
				for j = 1:m.spatDim(r)
					xx = repmat(xy(:,j),1,np);
					dx = xx-xx';
					di = di+dx.^2;
				end
				distance = sqrt(di);
				
				iWg = nan(np,np,alphagN);
				iWgR = nan(np,np,alphagN);
				detWg = nan(alphagN,1);
			case 'PGP'
				dim = m.spatDim(r);
				if isempty(m.xyKnot{r})
					if ~isnan(m.knotNumber(r))
						knotN = m.knotNumber(r);
					else
						error('HMSC: either knot locations or knot number must be specified for each spatial level with Gaussian predictive process');
					end
					xyMin = min(m.xy{r});
					intLim = max(m.xy{r})-xyMin;
					if dim==1
						hStep = ((1:R).^0)';
						hStep = repmat(hStep./sum(hStep),1,dim) .* repmat(intLim,R,1) ;
						hG0 = (cumsum(hStep)-hStep/2);
						hG = hG0';
						nK = size(hG,1);
					else
						side = ((prod(intLim) / knotN))^(1/dim);
						R = ceil(intLim / side);
						hStep = intLim ./ R;
						des = fullfact(R);
						hG0 = cell(1,dim);
						nK = prod(R);
						hG = nan(nK, dim);
						for d=1:dim
							hG0{d} = hStep(d)*((1:R(d))-0.5);
							hG(:,d) = hG0{d}(des(:,d));
						end
					end
					xyK = repmat(xyMin,nK,1) + hG;
					m.xyKnot{r} = xyK;
				else
					xyK = m.xyKnot{r};
					nK = size(xyK,1);
				end
				
				xy = m.xy{r};
				di = zeros(np,nK); % calculate the spatial distances
				for j = 1:dim
					xx1 = repmat(xy(:,j),1,nK);
					xx2 = repmat(xyK(:,j),1,np);
					dx = xx1-xx2';
					di = di+dx.^2;
				end
				di = sqrt(di);
				di12 = di;

				di = zeros(nK); % calculate the spatial distances
				for j = 1:dim
					xx = repmat(xyK(:,j),1,nK);
					dx = xx-xx';
					di = di+dx.^2;
				end
				di = sqrt(di);
				di22 = di;
				
				idDg = nan(np,alphagN);
				idDW12g = nan(np,nK,alphagN);
				Fg = nan(nK,nK,alphagN);
				iFg = nan(nK,nK,alphagN);
				detDg = nan(alphagN,1);
			case 'NNGP'
				iWg = cell(1,alphagN);
				iWgR = cell(1,alphagN);
				detWg = nan(alphagN,1);
				indNN = knnsearch(xy,xy,'K',m.nnNumber(r)+1);
				indNN = sort(indNN(:,2:end),2);
				indices = cell(2,np-1);
				distCell = cell(1,np);
				for i = 2:np
					ind = indNN(i,:);
					ind = ind(ind<i);
					if ~isempty(ind)
						indices{1,i} = i*ones([1,length(ind)]);
						indices{2,i} = ind;
						distCell{i} = squareform(pdist(xy([ind,i],:)));
					end
				end
		end
		
		
		for ag=1:alphagN
			alpha = alphapw(ag,1);
			switch m.spatialMethod{r}
				case 'GP'
					if alpha < 1e-5
						W = eye(np);
					else
						W = exp(-distance/alpha);
					end
					iW = invChol_mex(W);
					iWg(:,:,ag) = iW;
					iWgR(:,:,ag) = chol(iW);
					cholW = chol(W);
					detW = 2*sum(log(diag(cholW)));
					detWg(ag) = detW;
				case 'PGP'
					if alpha < 1e-5
						W22 = eye(nK);
						W12 = zeros(np,nK);
					else
						W22 = exp(-di22/alpha);
						W12 = exp(-di12/alpha);
					end
					iW22 = invChol_mex(W22);
					D = W12*iW22*W12';
					dD = 1-diag(D);
					
					LiW22 = chol(iW22,'lower');
					idD = dD.^-1;
					
					tmp0 = repmat(idD,1,nK);
					idDW12 = tmp0.*W12; % iD1*W12;
					F = W22 + W12'*idDW12;
					iF = invChol_mex(F);
					tmp2 = W12*LiW22;
					DS = tmp2'*(tmp0.*tmp2) + eye(nK);
					LDS = chol(DS, 'lower');
					detD = sum(log(dD))+2*sum(log(diag(LDS)));
					idDg(:,ag) = idD;
					idDW12g(:,:,ag) = idDW12;
					Fg(:,:,ag) = F;
					iFg(:,:,ag) = iF;
					detDg(ag) = detD;
				case 'NNGP'
					if alpha < 1e-5
						iW = speye(np);
						iWR = speye(np);
						detW = 1;
					else
						D = zeros(1,np);
						values = cell(1,np-1);
						D(1) = 1;
						for i = 2:np
							if ~isempty(indices{2,i})
								Kp = exp(-distCell{i}/alpha);
								values{i} = (Kp(1:end-1,1:end-1)\Kp(1:end-1,end))';
								D(i) = Kp(end,end) - Kp(end,1:end-1)*values{i}';
							else
								D(i) = 1;
							end
						end
						A = sparse(cell2mat(indices(1,:)),cell2mat(indices(2,:)),cell2mat(values),np,np);
						B = speye(np)-A;
						iWR = sparse(1:np,1:np,D.^-0.5)*B;
						iW = iWR'*iWR;
						detW = sum(log(D));
					end
					iWg{ag} = iW;
					iWgR{ag} = iWR;
					detWg(ag) = detW;
% 					if ag == alphagN
% 						[LiU,p,ord] = chol(iW,'lower');
% 						aaa
% 					end
					
			end
		end
		
		iWgA{r} = iWg;
		iWgRA{r} = iWgR;
		detWgA{r} = detWg;
		idDgA{r} = idDg;
		idDW12gA{r} = idDW12g;
		FgA{r} = Fg;
		iFgA{r} = iFg;
		detDgA{r} = detDg;
	end
end

end