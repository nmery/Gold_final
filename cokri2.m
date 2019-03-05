function [x0s,s,id,l,k,k0,dif,nbimag]=cokri2(x,x0,id,model,c,sv,itype,avg,ng)
%
% COKRI2 : This function is called from COKRI. The description for input and
%          output is given in COKRI. The only new variables are 'k0' which is
%          the right member matrix of the cokriging system and 'ng' which is
%          the total number of points for block discretization.
% Author: D. Marcotte
% Version 2.0  7/aug/14
% External subroutines: trans, means

x0s=[];s=[];l=[];k=[];k0=[];nc=0;dif=[];nbimag=[];

% definition of some constants

[n,t]=size(x);
[rp,p]=size(c);
r=rp/p;
[m,d]=size(x0);
cx=[x(:,1:d);x0];

% if no samples found in the search radius, return NaN
if n==0
  x0s=NaN*ones(m/ng,p);
  s=NaN*ones(m/ng,p);
  return
end


% calculation of left covariance matrix K and right covariance matrix K0

k=covardm(x(:,1:d),x(:,1:d),model,c);
k0=covardm(x(:,1:d),x0,model,c);
k00=covardm(x0,x0,model,c);

% constraints are added according to cokriging type

if itype==99,

  % the system does not have to be solved
  return
end


if itype==1.5;
  
  % krigeage oridinaire avec 1 seule contrainte sur la var. auxiliaire
  'arrêt'
end

%nbimag=0;
  
if itype==1.7;         % Cressies estimator (adapted here only for the kriging case, one point at a time, no missing value, no cokriging)
    
    ki=inv(k);         % inverse of variance-covariance matrix
    ski=sum(ki(:));    % sum of all elements in ki
    k0c=mean(k0,2);     % compute point-block covarianc  
    l=ki*k0c;
    k0l=k0c'*l;
    mu2=sqrt(l'*k0c*ski-sum(l).^2)/sqrt(ski*sv-1);   % Eq. 3.13
    mu1=-mu2/ski+sum(l)/ski;                   % Eq. 3.12
    lc=[1/mu2*k0c'*ki-mu1/mu2*sum(ki)];      % Eq. 3.11  
    l=lc';
   % dif=l'*k*l-sv;

    if sum(abs(imag(l)))>0 
   % k2=[k k0;k0' k00];
   % [m,nnn]=size(k);
    kg=[k ones(n,1);ones(1,n) 0]; % left hand kriging matrix
    kd=[mean(k0,2);1];    % right hand kriging vector
    l1=kg\kd;                % kriging weights
    lmm=l1(1:n);              % kriging weights without Lagrange multiplier
    l=lmm;
    nbimag=1;
    else
        nbimag=0;
    end
  %  x0s=l(2:end)'*x(2:end,end);
  %  s=sv+lc(2:end)'*k(2:end,2:end)*lc(2:end)-2*lc(2:end)'*k0c(2:end);
    
    x0s=l'*x(:,end);
    s=sv+l'*k*l-2*l'*k0c;
    l=[l;mu1;mu2];
  %  dif=l'*k*l-sv;
         
   return
 
end
    
if itype==2,
  
  % cokriging with one non-bias condition (Isaaks and Srivastava, 1990, p.410)
  
  k=[k,ones(n*p,1);ones(1,n*p),0]; k0=[k0;ones(1,m*p)]; nc=1;
elseif itype>=3,
  
  % ordinary cokriging (Myers, Math. Geol, 1982)
  
  t=kron(ones(1,n),eye(p));
  k=[k,t';t,zeros(p,p)];
  k0=[k0;kron(ones(1,m),eye(p))]; nc=p;
  nbimag=0;
  
  %   % cokriging with one non-bias condition in the z direction
  
  if itype == 3.5,
    
    t=kron(cx(1:n,d),eye(p));
    k=[k,[t;zeros(p,p)];[t',zeros(p,p+p)]];
    t=kron(cx(n+1:n+m,d)',eye(p));
    k0=[k0;t];
    nc=nc+p;
  end;
  
   
  if itype >=4,
    
    % universal cokriging ; linear drift constraints
    nca=p*d;
    t=kron(cx(1:n,:),eye(p));
    k=[k,[t;zeros(p,nca)];[t',zeros(nca,nc+nca)]];
    t=kron(cx(n+1:n+m,:)',eye(p));
    k0=[k0;t];
    nc=nc+nca;
  end;
  if itype==5,
    
    % universal cokriging ; quadratic drift constraints
    
    nca=p*d*(d+1)/2;
    cx2=[];
    for i=1:d,
      for j=i:d,
        cx2=[cx2,[cx(:,i).*cx(:,j)]];
      end
    end
    t=kron(cx2(1:n,:),eye(p));
    k=[k,[t;zeros(nc,nca)];[t',zeros(nca,nc+nca)]];
    t=kron(cx2(n+1:n+m,:)',eye(p));
    k0=[k0;t];
    nc=nc+nca;
  end
end

% columns of k0 are summed up (if necessary) for block cokriging

m=m/ng;
t=[];
for i=1:m,
  for ip=1:p,
    j=ng*p*(i-1)+ip;
    t=[t,means(k0(:,j:p:i*ng*p)')'];
  end
end
k0=t;

t=x(:,d+1:d+p);
if itype<3,
  
  % if simple cokriging or cokriging with one non bias condition, the means
  % are substracted
  
  t=(t-ones(n,1)*avg)';
else
  t=t';
end

% removal of lines and columns in K and K0 corresponding to missing values

z=zeros(n*p,1);
z(:)=t;
iz=~isnan(z);
iz2=[iz;ones(nc,1)]~=0;
nz=sum(iz);

% if no samples left, return NaN

if nz==0,
  x0s=nan;
  s=nan;
  return;
else
  k=k(iz2,iz2');
  k0=k0(iz2,:);
  id=id(iz,:);
  
  % check there is at least one sample for each variable
  
  [n1,m1]=size(k);
  for j=n1:-1:n1-nc+1;
    if k(j,:)==zeros(1,m1)
      k=k(1:j-1,1:j-1);
      k0=k0(1:j-1,:);
      m1=m1-1;
    end
  end
  
  
  % solution of the cokriging system by gauss elimination
  
  l=k\k0;
  
  % calculation of cokriging estimates
  
  t2=l(1:nz,:)'*z(iz);
  t=zeros(p,m);
  t(:)=t2;
  
  % if simple or cokriging with one constraint, means are added back
  
  if itype<3,
    t=t'+ones(m,1)*avg;
  else
    t=t';
  end
  x0s=t;
  
  % calculation of cokriging variances
  
  s=kron(ones(m,1),sv);
  t=zeros(p,m);
  t(:)=diag(l'*k0);
  s=s-t';
end

function [k]=Covardm(x,x0,model,c)
% [k]=Covar(x,x0,model,c)
% Fonction pour calculer les covariances avec des mod`eles spécifiées comme avec cokri
% La fonction calcule pour toutes les positions de x (n1) avec toutes les positions de x0 (n2)  K est donc n1xn2
% auteur D. Marcotte avril 2002

% here we define the equations for the various covariograms. Any new model
% can be added here.
k=[];
Gam=['h==0                                                               '; %nugget
  'exp(-h)                                                            '; %exponential
  'exp(-(h).^2)                                                       '; %gaussian
  '1-(1.5*min(h,1)/1-.5*(min(h,1)/1).^3)                              '; %spherical
  '1-h                                                                '; %linear
  '1-(7*min(h,1).^2-8.75*min(h,1).^3+3.5*min(h,1).^5-0.75*min(h,1).^7)'; %modele cubique
  '(h.^2).*log(max(h,eps))                                            '; %spline plaque mince
  '(h.^2+1).^(-0.5)                                                   '; %modèle gravimétrique (Cauchy avec b=0.5)
  '(h.^2+1).^(-1.5)                                                   '; %modele magnétique (Cauchy avec b=1.5)
  'sin(max(eps,h*2*pi))./max(eps,h*2*pi)                              '; %effet de trou sinusoidal
  'cos(h*2*pi)                                                        '; %effet de trou cosinusoidal
  '1-(1.5*min(h,1)/1-.5*(min(h,1)/1).^3)+1-h                          '];%spherique+lineaire


% definition of some constants


[n1,d]=size(x); % d dimension de l'espace
[n2,d]=size(x0);
[rp,p]=size(c);
r=rp/p;  % nombre de structures
cx=[x(:,1:d);x0];
nm=size(model,2);

% ne pas permettre des portées de 0 en input pour éviter divisions par 0
if nm>2
  model(:,2:1+d)=max(model(:,2:1+d),100*eps);
else
  model(:,2)=max(model(:,2),100*eps);
end
% calculer les covariances

k=zeros(n1*p,n2*p);
for i=1:r,
  
  % calculation of matrix of reduced rotated distances H
  
  [t1]=trans(x(:,1:d),model,i);
  [t2]=trans(x0,model,i);
  h=0;
  for id=1:d
    h=h+(t1(:,id)*ones(1,n2)-ones(n1,1)*t2(:,id)').^2;
  end
  h=sqrt(h);
  ji=(i-1)*p+1; js=i*p ;
  
  % evaluation of the current basic structure
  g=eval(Gam(model(i,1),:));
  k=k+kron(g,c(ji:js,:));
end

