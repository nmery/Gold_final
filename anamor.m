function [z]=anamor(y)
% fonction pour transformer des données en valeur gaussienne
% Syntaxe [z]=anamor(y)
% syntaxe : [z]=anamor(y)
% y est une matrice n x nc+1, nc est le nombre de coordonnées (nc>= 1)
%
[n,p]=size(y);
nc=p-1;
n2=length(unique(y(:,end)));

if n2<n  % there are one or more equal values. We use a kriging to break the equalities
    d=max(y(:,1:nc))-min(y(:,1:nc));
    dmax=max(d);
    
    %ys=cokri(y,zeros(1,nc),[1 1; 5 1],[0.5; 1/dmax],3,0,ones(1,nc),ones(1,nc),1,5,dmax/5,1);
    ys=y+rand(length(y),1)*0.0000000000001;
else
    ys=zeros(n,1);
end

[yt,id]=sortrows([y(:,end) ys(:,end)],[1 2]);

[t,id2]=sort(id);

rang=[1:n]';

z=norminv(rang/(n+1),0,1);

z=[y z(id2)];

