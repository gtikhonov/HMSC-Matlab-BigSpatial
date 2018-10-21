function betav = update_betav(X,Xr,Xs,Xv,z,beta,nf,pi,ncr,eta,lambda,etas,lambdas,qv,mu0v,V0v,sigma,outlierSpecies,spatial,factorcov,speciesX,includeXs)
[ny ns]=size(z);
ncv = size(Xv,2);
if speciesX
	nc=size(X{1},2);
else
	nc=size(X,2);
end
if length(spatial)==0
	np = [];
	nr = 0;
else
	[ny nr]=size(pi);
	np = max(pi);
end
S = z;
if speciesX
	for i=1:ns
		S(:,i) = S(:,i) - X{i}*beta(:,i);
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
for r=1:nr
	eta1=eta{r};
	lambda1=lambda{r};
	if factorcov(r) > 0
		Xr1 = Xr{r};
		for k=1:ncr(r)
			if factorcov(r) == 1
				XrEta = repmat(Xr1(:,k), 1, nf(r)).*eta1;
				S = S-XrEta(pi(:,r),:)*lambda1(:,:,k);
			else
				for j=1:ns
					XrEta = repmat(Xr1{j}(:,k), 1, nf(r)).*eta1;
					S(:,j) = S(:,j) - XrEta(pi(:,r),:)*lambda1(:,j,k);
				end
			end
		end
	else
		S = S-eta1(pi(:,r),:)*lambda1;
	end
end

betav =  zeros(ncv,ns);
for j=1:ns
    X1 = Xv.*repmat(qv(:,j)',ny,1);
    iVv = inv(V0v)+(1/sigma(j,j))*X1'*X1;
    Vv = inv(iVv);
    Vv = (Vv+Vv')/2;
    mbev = Vv*(inv(V0v)*mu0v+(1/sigma(j,j))*X1'*S(:,j));
    betav(:,j)=mvnrnd(mbev,Vv)';
end

