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
     '1-h.^2./(1+h.^2)                                                   '];%christakos 1984


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
 