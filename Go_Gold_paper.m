% Load datasets
load('ddh.dat');        % ddh comp=3, cap=30
load('real_tb.dat');    % rc simulated by TB
load('BM_uc.dat');      % Panel model for UC
load('ddh_dec.dat');    % declustered ddh comp=3, cap=30
load('ddh_rc.dat');     % ddh dataset based on rc (by ranking)

% Define datasets
x=ddh;
x_dec=ddh_dec;
G=real_tb(:,1:3);
R=real_tb;
G_UC=BM_uc;

% Parameters
block=[5 5 5];
nd=[3 3 3];
ival=0;
nk_ok=1;
nk_ck=3;
nk_uc=50;
nd_uc=[5 5 5];
rad=99999;
ntok=1;
avg=mean(x(:,4));
smu=[5 5 5];
panneau=[50 50 30];
vc=[0 0.05:0.01:1.5];
dec=1;      % dec=0 non-declus histogram for UC / dec=1 declus histogram for UC

%cas=1 best case ddh comp=3
%cas=2 best case using rc as ddh samples

vcas=[1 2];

for i=1:length(vcas)
    cas=vcas(i)
    switch cas
        case 1
            model=[1 1;4 25;4 60];
            c=[ 1.2;0.35;0.05];
            
        case 2
            x=ddh_rc;
            dec=0;
            nk_ck=1;
            model=[1 1;4 25;4 60];
            c=[ 1.2;0.35;0.05];
            
    end
    [stat,x0s_ok,x0s_ck,ton_uc,mean_uc]=okckuc(x,x_dec,G,R,model,c,block,nd,nd_uc,ival,nk_ok,nk_ck,nk_uc,rad,ntok,avg,smu,panneau,vc,G_UC,cas,dec);
    v_sum(i,:)=stat.summary;
    stats_real=[mean(real_tb(:,4)),var(real_tb(:,4)),std(real_tb(:,4))]
    stats_ok=stat.OK(2:end)
    stats_ck=stat.CK(2:end)
    stats_uc=mean_uc(1,1)
end