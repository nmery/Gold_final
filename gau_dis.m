
function zv=gau_dis(r,tab1);
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