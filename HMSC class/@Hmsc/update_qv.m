function qv = update_qv(X,Xr,Xs,Xv,z,beta,nf,pi,ncr,eta,lambda,etas,lambdas,betav,qv,q0v,sigma,outlierSpecies,spatial,factorcov,speciesX,includeXs)
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

for j=1:ns
    X1 = Xv.*repmat(qv(:,j)',ny,1);
    eps = S(:,j)-X1*betav(:,j);
    li = -(eps'*eps)/(2*sigma(j,j));
    for i=1:ncv
        if ~ismember(q0v(i,j),[0 1])
            mult=-betav(i,j);
            w=q0v(i,j);
            if qv(i,j)==0
                mult=-mult;
                w=1-w;
            end
            nw=1-w;
            dL = mult*Xv(:,i);
            neps = eps-dL;
            nli = -(neps'*neps)/(2*sigma(j,j));
            mini=min(li,nli);
            eli=exp(li-mini);
            enli=exp(nli-mini);
            if rand<(nw*enli)/(nw*enli+w*eli)
                qv(i,j)=1-qv(i,j);
                li=nli;
            end
        end
    end
end
