function sampleMCMC(m, nRuns, append, startPar, verb, calculateDataParametersFlag, anVec)

speciesX = m.speciesX;
phylogeny = m.phylogeny;
includeXs = m.includeXs;
includeXv = m.includeXv;
outlierSpecies = m.outlierSpecies;
spatial = m.spatial;
spatialMethod = m.spatialMethod;
factorCov = m.factorCov;

f0 = m.f0;
V0 = m.V0;
Ugamma = m.Ugamma;
iUgamma = invChol_mex(Ugamma);
mgamma = m.mgamma;
asigma = m.asigma;
bsigma = m.bsigma;

mu0v=m.mu0v;
V0v=m.V0v;
q0v=m.q0v;

nur = m.nur;
a1r = m.a1r;
b1r = m.b1r;
a2r = m.a2r;
b2r = m.b2r;
nus = m.nus;
a1s = m.a1s;
b1s = m.b1s;
a2s = m.a2s;
b2s = m.b2s;
nu = m.nu;

% nc = m.nc;
ns = m.ns;
% nt = m.nt;
ncr = m.ncr;
np = m.np;
dist = m.dist;

Xr = m.Xr;
Xs = m.Xs;
Xv = m.Xv;
T = m.T;
Y = m.Y;
pi = m.pi;
rhopw = m.rhopw;
alphapw = m.alphapw;

if calculateDataParametersFlag
   [iQg,detQg,iWg,iWgR,detWg,idDg,idDW12g,Fg,iFg,detDg] = m.calculateDataParameters();
   % [iQg,detQg,iWg,iWgR,detWg,Cg,taug,iICCg] = m.calculateDataParameters();
   m.iQg = iQg;
   m.detQg = detQg;
   m.iWg = iWg;
   m.iWgR = iWgR;
   m.detWg = detWg;
   m.idDg = idDg;
   m.idDW12g = idDW12g;
   m.Fg = Fg;
   m.iFg = iFg;
   m.detDg = detDg;
else
   iQg = m.iQg;
   detQg = m.detQg;
   iWg = m.iWg;
   iWgR = m.iWgR;
   detWg = m.detWg;
   idDg = m.idDg;
   idDW12g = m.idDW12g;
   Fg = m.Fg;
   iFg = m.iFg;
   detDg = m.detDg;
end

if append == false
	m.repPar = [];
	m.repN = 0;
	startRunN = 0;
	[X,covScale,Xs,sCovScale,Xv,vCovScale,Xr,factorCovScale, beta,gamma,iV,sigma,rho,ph,z,nf,lambda,eta,delta,psijh,alpha,nfs,lambdas,etas,deltas,psijhs,betav,qv] = m.computeInitialValues(startPar);
	% clear the directory
else
	startRunN = length(m.repPar);
	flag = 0;
	if isempty( m.repPar{m.repN} )
		flag = 1;
		m.postFileToRam(m.repN);
	end
	lastPar = m.repPar{m.repN}{end};
	if flag
		m.postRamClear();
	end
	[X,covScale,Xs,sCovScale,sCovScale,Xv,Xr,factorCovScale,beta,gamma,iV,sigma,rho,ph,z,nf,lambda,eta,delta,psijh,alpha,nfs,lambdas,etas,deltas,psijhs,betav,qv] = m.computeInitialValues(lastPar);
end
tic;
if verb > 0
	fprintf('Sampling starts.\n');
end
for run = startRunN+(1:nRuns)
	if verb > 1
		fprintf('SAVE REPLICATE: %d (%d/%d) --- %d\n',run,run-startRunN,nRuns,round(toc));
	end
	parVec = cell(1, length(m.samVec) );
	parK = 1;
	an = anVec(run);
	for mcmc = 1:m.samVec(end)
      z = Hmsc.update_z(z,X,Xr,Xs,Xv,beta,eta,lambda,etas,lambdas,betav,qv,Y,nf,pi,ncr,dist,sigma,speciesX,includeXs,includeXv,factorCov,an);
      eta = Hmsc.update_eta(z,X,Xr,Xs,Xv,beta,sigma,eta,alpha,lambda,etas,lambdas,betav,qv,nf,pi,ncr,spatial,spatialMethod,factorCov,iWg,idDg,idDW12g,Fg,speciesX,includeXs,includeXv,an);
      [beta,lambda] = Hmsc.update_betaLambda(Y,X,Xr,Xs,Xv,z,beta,nf,pi,ncr,eta,lambda,psijh,delta,etas,lambdas,betav,qv,...
          gamma,sigma,T,ph,iV,rho,phylogeny,iQg,detQg,outlierSpecies,spatial,factorCov,speciesX,includeXs,includeXv,an);
%       if run<=1
%           for r = 1:m.nr
%               etaMean = repmat(mean(eta{r},1),m.np(r),1);
%               etaSd = repmat(std(eta{r},1),m.np(r),1);
%               eta{r} = (eta{r}-etaMean)./repmat(etaSd,m.np(r),1);
%               lambda{r} = lambda{r}.*(etaSd');
%               beta(1,:) = beta(1,:) + etaMean*lambda{r};
%           end
%       end
      alpha = Hmsc.update_alpha(eta,alpha,nf,pi,spatial,spatialMethod,iWg,iWgR,detWg,idDg,idDW12g,iFg,detDg,alphapw,an);
      sigma = Hmsc.update_sigma(X,Xr,Xs,Xv,sigma,beta,eta,lambda,etas,lambdas,betav,qv,nf,pi,ncr,dist,asigma,bsigma,z,spatial,factorCov,speciesX,includeXs,includeXv,an);
      [gamma,iV] = Hmsc.update_gamma_V(T,beta,gamma,ph,rho,V0,f0,iUgamma,mgamma,phylogeny,iQg,outlierSpecies,anVec(run));
      [psijh,delta] = Hmsc.update_lambda_priors(nf,nur,a1r,a2r,b1r,b2r,psijh,delta,lambda,factorCov,anVec(run));
      
      if m.nr > 0
          if (run<m.adapt(1) || (run==m.adapt(1) && mcmc<=m.adapt(2)) ) && isequal(m.fixNf, false)
              [nf,lambda,eta,alpha,psijh,delta] = Hmsc.update_nf(nf,ncr,ns,mcmc,np,nur,a2r,b2r,lambda,eta,alpha,psijh,delta,spatial);
          end
      end
      
      if any(mcmc==m.samVec)
          if verb > 2
              fprintf(strcat('save replicate:', int2str(run), ',run:', int2str(run-startRunN), ',iteration round:',int2str(mcmc),'---', num2str(round(toc)),'\n'));
          end
          p = m.addToOutput(beta,gamma,sigma,rho,iV,lambda,eta,alpha,ph,nf,lambdas,etas,nfs,betav,qv,psijh,delta,psijhs,deltas,covScale,sCovScale,vCovScale,factorCovScale);
          parVec{parK} = p;
          parK = parK + 1;
		end
	end
	if m.ramPost == true
		m.repPar{run} = parVec;
	else
		m.repPar{run} = [];
	end
	if m.stfPost == true
		m.stfMCMCRun(parVec, run);
	end
	m.repN = m.repN + 1;
end
if verb>0
	toc;
end
end