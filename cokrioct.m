function [x0s,s,sv,idout,l,k,k0,nbimag2]=cokrioct(x,x0,model,c,itype,avg,block,nd,ival,nk,rad,ntok)
%[x0s,s,sv,idout,l,k,k0]=oct1(x,x0,model,c,itype,avg,block,nd,ival,nk,rad,ntok%)
%
% Pré-processeur de cokri pour la recherche par octant
% Ici, nk est le nombre de points par côté, quadrant ou octant (en 1D,2D ou 3D)
% En 4D ou plus, le programme fait une recherche circulaire simple.
% ntok est toujours égl à 1 lors de l'appel à cokridir. cokridir est une version modifiée
% de cokri ou l'on a enlevé la recherche de points.
%
% 
%
% COKRI performs point or block cokriging in D dimensions (any integer)
%       of P variables (any integer) with a combination of K basic models
%       (any integer).
%
% Syntax:
% [x0s,s,sv,id,l]=cokri(x,x0,model,c,itype,avg,block,nd,ival,nk,rad,ntok)
%
% Input description:
%   x: The n x (p+d) data matrix. This data matrix can be imported from an
%      existing ascii file. Missing values are coded 'nan' (not-a-number).
%   x0: The m x d matrix of coordinates of points to estimate.
%   model: Each row of this matrix describes a different elementary structure.
%          The first column is a code for the model type, the d following
%          columns give the ranges along the different coordinates and the
%          subsequent columns give rotation angles (a maximum of three).
%	   For more details on how to specify rotations, type help trans.
%          The codes for the current models are:
%             1: nugget effect
%             2: exponential model
%             3: gaussian model
%             4: spherical model
%             5: linear model
%          Note: a linear model is specified by arbitrary ranges and a sill
%                such that sill/range gives the desired slope in the direction
%                condidered.
%   c: The (rp x p) coefficient matrix of the coregionalization model.
%      Position (i,j) in each submatrix of size p x p give the sill of the
%      elementary component for each cross-variogram (variogram) between
%      variable i and variable j.
%   itype: Code to indicate which type of cokriging is to be performed:
%             1:  simple cokriging
%             2:  ordinary cokriging with one nonbias condition
%                 (Isaaks and Srivastava).
%             3:  ordinary cokriging with p nonbias condition.
%             4:  universal cokriging with drift of order 1.
%             5:  universal cokriging with drift of order 2.
%             99: cokriging is not performed, only sv is computed.
%   block: Vector (1 x d), giving the size of the block to estimate;
%          any values when point cokriging is required.
%   nd: Vector (1 x d), giving the discretization grid for block cokriging;
%       put every element equal to 1 for point cokriging.
%   ival: Code for cross-validation.
%             0:  no cross-validation
%             1:  cross-validation is performed  by removing one variable at a
%                 time at a given location.
%             2:  cross-validation is performed by removing all variables at a
%                 given location.
%   nk: Number of nearest neighbors in x matrix to use in the cokriging
%       (this includes locations with missing values even if all variables
%       are missing).
%   rad: Search radius for neighbors.
%   ntok: Points in x0 will be kriged by groups of ntok grid points.
%         When ntok>1, the search will find the nk nearest samples within
%         distance rad from the current ntok grid points centroid.
%????????????????????????????????????????????????????????????????????????????????????????????????????????????
% Output description:
%
%   For the usual application, only x0s and s are required and the other
%   output matrices may be omitted.
%
%   x0s: m x (d+p) matrix of the m points (blocks) to estimate by the
%        d coordinates and p cokriged estimates.
%   s: m x (d+p) matrix of the m points (blocks) to estimate by the
%      d coordinates and the p cokriging variances.
%   sv: 1 x p vector of variances of points (blocks) in the universe.
%   id: (nk x p) x 2 matrix giving the identifiers of the lambda weights for
%       the last cokriging system solved.
%   l: ((nk x p) + nc) x (ntok x p) matrix with lambda weights and
%      Lagrange multipliers of the last cokriging system solved.
%
% Syntax:
% [x0s,s,sv,id,l]=cokri(x,x0,model,c,itype,avg,block,nd,ival,nk,rad,ntok)
%
% Author: D. Marcotte
% Version 2.0  97/aug/14 (matlab4 and 5 compatible)
% External subroutines: cokri2, trans, means

[nr,ncol]=size(x);
[npts,bb]=size(x0);
p=size(c,2);
nvar=ncol-bb;
x0s=[x0 nan*ones(npts,p)];
s=[x0 nan*ones(npts,p)];
ntok=1;
rad2=rad^2;
nbimag2=[zeros(npts,p)];


% Cas 1D

if bb==1;
  for i=1:npts,  % on boucle sur chaque point à estimer
    d=abs(x0(i,1)-x(:,1));
    id1=d<rad & x0(i,1)<x(:,1);xid1=x(id1,:);
    id2=d<rad & x0(i,1)>=x(:,1);xid2=x(id2,:); 
    [q1,t1]=sort(d(id1,:));n1=sum(id1);
    [q2,t2]=sort(d(id2,:));n2=sum(id2);
    xc=[xid1(t1(1:min(n1,nk)),:);xid2(t2(1:min(n2,nk)),:)];
    nt=size(xc,1);
    plot(x0(i,1),x0(i,1),'*r',x(:,1),x(:,1),'+k',xc(:,1),xc(:,1),'ob')
    if nt>0
      [x0s(i,:),s(i,:),sv,idout,l,k,k0,nbimag2]=cokridir(xc,x0(i,:),model,c,itype,avg,block,nd,ival,nt,rad,ntok);
    else
      x0s(i,:)=nan;
      s(i,:)=nan;
    end
  end
  
  % cas 2D
  
elseif bb==2;
  if ival>=1
    npts=size(x,1);
  end
  
    for ii=1:npts, % on boucle sur les points à estimer
      if ival>=1;
    
      x0=x(ii,1:2);
      i=1;
    else
     i=ii;  
    end
    dt=(x0(i,1)-x(:,1)).^2+(x0(i,2)-x(:,2)).^2;
    xt=x;
    id=~isnan(x(:,bb+1));
    x=xt(id,:);
    d=dt(id);
    id1=d<rad2 & x0(i,1)<x(:,1) & x0(i,2)<=x(:,2);xid1=x(id1,:);
    id2=d<rad2 & x0(i,1)>=x(:,1) & x0(i,2)<x(:,2);xid2=x(id2,:);
    id3=d<rad2 & x0(i,1)>x(:,1) & x0(i,2)>=x(:,2);xid3=x(id3,:);
    id4=d<rad2 & x0(i,1)<=x(:,1) & x0(i,2)>x(:,2);xid4=x(id4,:);
    [q1,t1]=sort(d(id1,:));n1=sum(id1);
    [q2,t2]=sort(d(id2,:));n2=sum(id2);
    [q3,t3]=sort(d(id3,:));n3=sum(id3);
    [q4,t4]=sort(d(id4,:));n4=sum(id4);
    xc=[xid1(t1(1:min(n1,nk)),:);xid2(t2(1:min(n2,nk)),:);...
        xid3(t3(1:min(n3,nk)),:);xid4(t4(1:min(n4,nk)),:)];
    
    if nvar>1;
      for ivar=2:nvar
        id=~isnan(x(:,bb+ivar));
        x=xt(id,:);
        d=dt(id,:);
        id1=d<rad2 & x0(i,1)<x(:,1) & x0(i,2)<=x(:,2);xid1=x(id1,:);
        id2=d<rad2 & x0(i,1)>=x(:,1) & x0(i,2)<x(:,2);xid2=x(id2,:);
        id3=d<rad2 & x0(i,1)>x(:,1) & x0(i,2)>=x(:,2);xid3=x(id3,:);
        id4=d<rad2 & x0(i,1)<=x(:,1) & x0(i,2)>x(:,2);xid4=x(id4,:);
        [q1,t1]=sort(d(id1,:));n1=sum(id1);
        [q2,t2]=sort(d(id2,:));n2=sum(id2);
        [q3,t3]=sort(d(id3,:));n3=sum(id3);
        [q4,t4]=sort(d(id4,:));n4=sum(id4);
        xc=[xc;xid1(t1(1:min(n1,nk)),:);xid2(t2(1:min(n2,nk)),:);...
            xid3(t3(1:min(n3,nk)),:);xid4(t4(1:min(n4,nk)),:)];
      end
      x=xt;
      if ival>=1;
        xc=[x(ii,:);xc];
      end
      
        % enlever les points apparaissant plus d'une fois
      nt=size(xc,1);jj=1;
      while jj<=nt;
        num=[jj+1:nt];
        id=xc(jj+1:end,1)==xc(jj,1) & xc(jj+1:end,2)==xc(jj,2);
        xc(num(id),:)=[];
        jj=jj+1;
        nt=nt-sum(id);
      end
      
    end
    
    nt=size(xc,1);
    if nt>0
   %         plot(x0(i,1),x0(i,2),'*r',x(:,1),x(:,2),'+k',xc(:,1),xc(:,2),'ob')
   %   disp(['Kriging point ',num2str(ii),' parmi ',num2str(npts)]);
        [x0s(i,:),s(i,:),sv,idout,l,k,k0,nbimag2]=cokridir(xc,x0(i,:),model,c,itype,avg,block,nd,ival,nt,rad,ntok);

    else
      x0s(i,:)=nan;
      s(i,:)=nan;
    end
  end
  % cas 3D      
elseif bb==3;
  
  for i=1:npts,  % Boucler sur les points à estimer
    
    d=(x0(i,1)-x(:,1)).^2+(x0(i,2)-x(:,2)).^2+(x0(i,3)-x(:,3)).^2;
    id1=d<rad2 & x0(i,1)<=x(:,1) & x0(i,2)<=x(:,2) & x0(i,3)<=x(:,3);xid1=x(id1,:);
    id2=d<rad2 & x0(i,1)>x(:,1) & x0(i,2)<=x(:,2) & x0(i,3)<=x(:,3);xid2=x(id2,:);
    id3=d<rad2 & x0(i,1)>x(:,1) & x0(i,2)>x(:,2) & x0(i,3)<=x(:,3);xid3=x(id3,:);
    id4=d<rad2 & x0(i,1)<=x(:,1) & x0(i,2)>x(:,2) & x0(i,3)<=x(:,3);xid4=x(id4,:);
    id5=d<rad2 & x0(i,1)<=x(:,1) & x0(i,2)<=x(:,2) & x0(i,3)>x(:,3);xid5=x(id5,:);
    id6=d<rad2 & x0(i,1)>x(:,1) & x0(i,2)<=x(:,2) & x0(i,3)>x(:,3);xid6=x(id6,:);
    id7=d<rad2 & x0(i,1)>x(:,1) & x0(i,2)>x(:,2) & x0(i,3)>x(:,3);xid7=x(id7,:);
    id8=d<rad2 & x0(i,1)<=x(:,1) & x0(i,2)>x(:,2) & x0(i,3)>x(:,3);xid8=x(id8,:);
    [q1,t1]=sort(d(id1,:));n1=sum(id1);
    [q2,t2]=sort(d(id2,:));n2=sum(id2);
    [q3,t3]=sort(d(id3,:));n3=sum(id3);
    [q4,t4]=sort(d(id4,:));n4=sum(id4);
    [q5,t5]=sort(d(id5,:));n5=sum(id5);
    [q6,t6]=sort(d(id6,:));n6=sum(id6);
    [q7,t7]=sort(d(id7,:));n7=sum(id7);
    [q8,t8]=sort(d(id8,:));n8=sum(id8);
    xc=[xid1(t1(1:min(n1,nk)),:);xid2(t2(1:min(n2,nk)),:);...
        xid3(t3(1:min(n3,nk)),:);xid4(t4(1:min(n4,nk)),:);...
        xid5(t5(1:min(n5,nk)),:);xid6(t6(1:min(n6,nk)),:);...
        xid7(t7(1:min(n7,nk)),:);xid8(t8(1:min(n8,nk)),:)];
    nt=size(xc,1);
%    plot3(x0(i,1),x0(i,2),x0(i,3),'*r',x(:,1),x(:,2),x(:,3),'+k',xc(:,1),xc(:,2),xc(:,3),'ob')
    
    if nt>0
      [x0s(i,:),s(i,:),sv,idout,l,k,k0,nbimag2(i)]=cokridir(xc,x0(i,:),model,c,itype,avg,block,nd,ival,nt,rad,ntok);
       % nbimag3=[nbimag3;nbimag2];
    else
      x0s(i,:)=nan;
      s(i,:)=nan;
    end
  end
end


