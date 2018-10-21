function [gamma,iV] = update_gamma_V(T,beta,gamma,ph,rho,V0,f0,iUgamma,mgamma,phylogeny,iQg,outlierspecies,an)
[ns nt]=size(T);
[nc ns]=size(beta);

if phylogeny
	iQ=iQg(:,:,rho);
else
	iQ=eye(ns);
end

Ts = T;
betas=beta;
if outlierspecies
	for i=1:ns
		Ts(i,:)=sqrt(ph(i))*T(i,:);
		betas(:,i)=sqrt(ph(i))*beta(:,i);
	end
end

E=betas'-(Ts*gamma);
A=E'*iQ*E;
Vn=invChol_mex(A+V0);
% Vn=inv(A+V0);
Vn=(Vn+Vn')/2;
iV = wishrnd(Vn*an,nc+1+(f0+ns-nc-1)/an);

tmp=Ts'*iQ;
Ugammas=invChol_mex(iUgamma+kron(iV,tmp*Ts));
% Ugammas=inv(inv(Ugamma)+kron(iV,tmp*Ts));

Ugammas=(Ugammas+Ugammas')/2;
res=tmp*betas'*iV;
mgammas=Ugammas*(iUgamma*mgamma+res(:));
% mgammas=Ugammas*(inv(Ugamma)*mgamma+res(:));
gamma=mvnrnd(mgammas,Ugammas*an)';
gamma=reshape(gamma,nt,nc);
