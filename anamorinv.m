function [z]=anamorinv(y,w)
% fonction pour transformer les valeurs gaussiennes en valeur "vraie" 
% syntaxe : [z]=anamorinv(y,w)
% y est une variable à n x nc+2 entrées: [corrdonnées,valeur,valeur gaussienne]
% w est les données gaussiennes simulées (ns x 1)
%
% output: z, ns x 2, [valeur gaussienne, valeur originale]

[n,p]=size(y);
nc=p-2;
yminv=0;                  % valeur minimale vraie (changer au besoin)
yming=-7;                 % valeur minimale gaussienne (changer au besoin)
ymaxv=max(y(:,p-1))*1.1;  %le 1.1 est un facteur multiplicateur du maximum (changer au besoin)
ymaxg=7;                  % valeur maximale gaussienne

[yt,id]=sortrows([y(:,p-1) y(:,p)],[2]);

y=[yminv yming; yt(:,1) yt(:,2); ymaxv ymaxg];

z=interp1(y(:,2),y(:,1),w,'linear');

z=[w z]; %si w a des coordonnées

