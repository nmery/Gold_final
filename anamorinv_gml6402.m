function [z]=anamorinv_gml6402(tab,y)
% fonction pour transformer les valeurs gaussiennes y N(0,1) en valeur z en se servant du tableau d'équivalence tab 
% syntaxe : [z]=anamorinv_gml6402(tab,y)
% tab est un tableau n x 2 entrées: [valeur gaussienne, valeur z correspondante]
% y est un vecteur gaussien ns x 1
% le minimum possible sur z est indiqué être 0. Le maximum possible est
% pris 110% de la valeur maximale, ces valeurs extrêmes étant associées à
% -7 et 7 dans l'espace gaussien.
% output: z est ns x 2, [y,z]

[nc,~,~]=size(tab);
zmin=0;
%ymin=-7;                 % associer des valeurs extrêmes pour éviter des extrapolations
%zmax=max(tab(:,2))*1.1;  %le 1.1 est un facteur multiplicateur du maximum du tableau
%ymax=7;

[yt,id]=sortrows(tab,[1]); % s'assurer que tab est trié par valeurs croissantes

%tab=[ymin zmin ; yt ; ymax zmax];
tab=yt;

z=interp1(tab(:,1),tab(:,2),y,'linear','extrap');

z=[y z]; 

