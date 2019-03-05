function [Tc,mc,R,rV,rv,yV_pan,yV_fit,tabp,tabs]=cond_unif_v2(zk,model,c,tab,panneau,smu,vc)
% Syntaxe:
% zk : valeurs krigées des panneaux
% model, c: modèles de covariance (comme cokri) de la variable Z
% tab: tableau nx2 (y,z) donnant l'anamorphose sous forme graphique y est
% N(0,1)
% panneau : taille du panneau (1 x d) d la dimension de l'espace
% smu: taille du smu (1 x d) 
% vc: vecteur avec les teneurs de coupure
%
%Tc proportion tonnage dans le panneau
% mc teneur du minerai dans le panneau

d=size(panneau,2);
if d==1
    
    nd=100;
elseif d==2
    nd=[11,11];
elseif d==3;
    nd=[7,7,7];
  % nd=[20,10,10];
end
    
[x0s,s,s2p]=cokri([ones(1,d),1],zeros(1,d),model,c,3,0,panneau,nd,0,999,9999,1);
[x0s,s,s2v]=cokri([ones(1,d),1],zeros(1,d),model,c,3,0,smu,nd,0,999,9999,1);

fr=@(r,s2v,tab1,yv,u)trouve_r(r,s2v,tab1);
rV=fminbnd(fr,0,1,[],s2p,tab);
rv=fminbnd(fr,0,1,[],s2v,tab);
R=rV/rv;

tabp=anamor_p(rV,tab);    % calculer l'anamorphose au niveau des panneaux
tabs=anamor_p(rv,tab);   % anamorphose smu
tabp=[-70 0;tabp;70 1.5*max(tab(:,2))];
tabs=[-70 0;tabs;70 1.5*max(tab(:,2))];
yV_pan=interp1(tabp(:,2),tabp(:,1),zk(:,end),'linear','extrap');

p=[0.0005:0.001:1]';
u=norminv(p,0,1);
fyv=@(yV,zk,tab,u,R)trouve_yv(yV,zk,tab,u,R);
yV_fit=zeros(size(zk));
yc=interp1(tabs(:,2),tabs(:,1),vc,'linear','extrap');

for i=1:length(zk)
    
%    yv=u*sqrt(1-R^2)+R*yV(i);
%    zv=interp1(tabs(:,1),tabs(:,2),yv,'linear','extrap');
   yV=fminbnd(fyv,-70,70,[],zk(i),tabs,u,R);
   yV_fit(i)=yV;
   [fit,zvbar,zv]=trouve_yv(yV,zk(i),tabs,u,R);
   if abs(zvbar-zk(i))>0.01*zk(i)
       disp('attention impossible de rencontrer la teneur du panneau')
   end
   for j=1:length(vc);
       Tc(i,j)=1-normcdf((yc(j)-R*yV_fit(i))/sqrt(1-R^2));
       f = @(x,R,yV,tabs) interp1(tabs(:,1),tabs(:,2),x,'linear','extrap').*normpdf(x,R*yV,sqrt(1-R^2));
       mc(i,j)= integral(@(x)f(x,R,yV,tabs),yc(j),Inf)/Tc(i,j);
   end
end
disp('')


function tab=anamor_p(r,tab1);
% fonction pour calculer la valeur Zv avec la corrélation r et
% l'anamorphose tab1 utilisant les valeurs normales de Yv et u
% Syntaxe: zv=gau_dis(r,tab1,yv,u)
% r: corrélation entre les gaussiennes, coefficient de changement de support
% tab1: tableau donnant l'anamorphose na x 2, 1ere colonne, Y, 2e Z correspondant
% zv : valeur de bloc correspondant à Yv n1 x 1


u=norminv([0.0005:0.001:1]',0,1);   % définir les valeurs N(0,1)
yv=u;
ny=length(u);

y=r*kron(yv,ones(ny,1))+sqrt(1-r^2)*kron(ones(ny,1),u);
z=anamorinv_gml6402(tab1,y);
zv=mean(reshape(z(:,end),ny,ny));   % calculer la teneur moyenne de bloc
zv=zv(:);
tab=[yv(:) zv(:)];
tab(:,2)=tab(:,2)+eps*100*[1:size(yv,1)]'; % éliminer tout risque d'égalité




