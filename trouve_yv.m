function [fit,zvbar,zv]=trouve_yv(yV,zk,tab,u,R)
% fonction pour trouver le r permettant d'avoir la variance de bloc s2v
% r : r actuel
% s2v: variance de bloc désirée
% tab1: tableau décrivant l'anamorphose

yv=u*sqrt(1-R^2)+R*yV;
zv=interp1(tab(:,1),tab(:,2),yv,'linear','extrap');
zvbar=mean(zv);
fit=(zk-zvbar).^2;

