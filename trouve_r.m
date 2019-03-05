
function [fit]=trouve_r(r,s2v,tab1)
% fonction pour trouver le r permettant d'avoir la variance de bloc s2v
% r : r actuel
% s2v: variance de bloc désirée
% tab1: tableau décrivant l'anamorphose

zv=gau_dis(r,tab1);
s2v_gau=var(zv);
fit=(s2v-s2v_gau).^2;