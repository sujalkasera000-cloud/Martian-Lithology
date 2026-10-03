function [kbr,mubr,vp,vs,ro2,k2] = berryscm_upper_halite_sc1(k,mu,asp,x,ro1,P_gas)
%BERRYSCM - Effective elastic moduli for multi-component composite
% using Berryman's Self-Consistent (Coherent Potential Approximation) method. 
%
%[KBR,MUBR,VP,VS,RO2,K2]=BERRYSCM(K,MU,ASP,X,RO1,P_GAS)
%	K,MU:       Bulk and shear moduli of the N constituent phases (K, MU, vectors of length N)
%	ASP:        Aspect ratio for the inclusions of the N phases < 1 for oblate spheroids; >1 for prolate spheroids.
%   X:          Fraction of each phase. Solid, then fluid phase. Sum(X) should be 1.
%   RO1:        Initial density of the dry rock.
%   P_GAS:      Fraction of gas in the pores as a decimal (e.g., 0.98 for 98% gas). 
%               The remaining pore space is left blank (vacuum).
%	KBR,MUBR:  	Effective bulk and shear moduli 
%
%   Original berryman function Written by T. Mukerji, 1997, then modified by Vashan
%

%% Calculate elastic moduli
k=k(:); mu=mu(:); asp=asp(:); x=x(:);
indx=find(asp==1); asp(indx)=0.99*ones(size(indx));
theta=zeros(size(asp)); fn=zeros(size(asp));  

obdx=find(asp<1);
theta(obdx)=(asp(obdx)./((1-asp(obdx).^2).^(3/2))).*...
             (acos(asp(obdx)) -asp(obdx).*sqrt(1-asp(obdx).^2));
fn(obdx)=(asp(obdx).^2./(1-asp(obdx).^2)).*(3.*theta(obdx) -2);

prdx=find(asp>1);
theta(prdx)=(asp(prdx)./((asp(prdx).^2-1).^(3/2))).*...
             (asp(prdx).*sqrt(asp(prdx).^2-1)-acosh(asp(prdx)));
fn(prdx)=(asp(prdx).^2./(asp(prdx).^2-1)).*(2-3.*theta(prdx));

ksc= sum(k.*x);
musc= sum(mu.*x);
knew= 0.;
tol=1e-6*k(1);
del=abs(ksc-knew);
niter=0;

while( (del > abs(tol)) && (niter<3000) ) 
	nusc=(3*ksc-2*musc)/(2*(3*ksc+musc));
	a=mu./musc -1; 
	b=(1/3)*(k./ksc -mu./musc); 
	r=(1-2*nusc)/(2*(1-nusc));
	
    f1=1+a.*((3/2).*(fn+theta)-r.*((3/2).*fn+(5/2).*theta-(4/3)));
	f2=1+a.*(1+(3/2).*(fn+theta)-(r/2).*(3.*fn+5.*theta))+b.*(3-4*r);
	f2=f2+(a/2).*(a+3.*b).*(3-4.*r).*(fn+theta-r.*(fn-theta+2.*theta.^2));
	f3=1+a.*(1-(fn+(3/2).*theta)+r.*(fn+theta));
	f4=1+(a./4).*(fn+3.*theta-r.*(fn-theta));
	f5=a.*(-fn+r.*(fn+theta-(4/3))) + b.*theta.*(3-4*r);
	f6=1+a.*(1+fn-r.*(fn+theta))+b.*(1-theta).*(3-4.*r);
	f7=2+(a./4).*(3.*fn+9.*theta-r.*(3.*fn+5.*theta)) + b.*theta.*(3-4.*r);
    f8=a.*(1-2.*r+(fn./2).*(r-1)+(theta./2).*(5.*r-3))+b.*(1-theta).*(3-4.*r);
	f9=a.*((r-1).*fn-r.*theta) + b.*theta.*(3-4.*r);
	
    p=3*f1./f2; 
	q=(2./f3) + (1./f4) +((f4.*f5 + f6.*f7 - f8.*f9)./(f2.*f4));
	
    p=p./3; 
	q=q./5; 
%------------------------------------------------------------------------
	knew= sum(x.*k.*p)/sum(x.*p);
	munew= sum(x.*mu.*q)/sum(x.*q);
	
	del=abs(ksc-knew);
	ksc=knew;
	musc=munew;
	niter=niter+1;
end		
kbr=ksc; mubr=musc;

%% Density, elastic moduli, and porosity parameters
rofl1 = 0.020;                  % density of gas
kfl1 = 0;                       % bulk modulus of gas

% The remaining pore space is left blank (density = 0, bulk modulus = 0)
rofl2 = P_gas * rofl1;          % effective density of fluid (gas + blank space)
kfl2 = P_gas * kfl1;            % effective bulk modulus of fluid (gas + blank space)

k0 = k(1);                      % bulk moduli of solid mineral phase
phi = x(2);                     % porosity of rock
k1 = kbr;                       % dry bulk modulus

%% Modify the seismic velocities to account for fluid filled rock
% Density after partial gas fill (remaining porosity acts as a vacuum)
ro2 = ro1 - phi.*rofl1 + phi.*rofl2;                       

% Fluid substitution using Gassmann equation
a = k1./(k0-k1) - kfl1./(phi.*(k0-kfl1)) + kfl2./(phi.*(k0-kfl2)); 
k2 = k0.*a ./ (1+a);                                      % bulk modulus
mu2 = mubr;                                               % shear modulus does not change

% Calculate Seismic Velocities
vp = sqrt((k2+(4/3)*mu2)./ro2);                           % vp after fluid substitution
vs = sqrt(mu2./ro2);                                      % vs after fluid substitution

end