function [stat,x0s_ok,x0s_ck,ton_uc,mean_uc]=okckuc(x,x_dec,G,R,model,c,block,nd,nd_uc,ival,nk_ok,nk_ck,nk_uc,rad,ntok,avg,smu,panneau,vc,G_UC,cas,dec)

[~,x2]=size(x);

if x2==3
    
    if dec==0
        %OK
        itype=3;
        [x0s_ok,~,~,~,~,~,~,nbimag_ok]=cokrioct(x,G,model,c,itype,avg,block,nd,ival,nk_ok,rad,ntok);
        nbimag_ok=sum(nbimag_ok);
        neg_ok=length(find(x0s_ok(:,3)<0));
        nan_ok=length(find(isnan(x0s_ok)));
        
        x0s_ok(x0s_ok<0)=0;
        mean_ok=nanmean(x0s_ok(:,3));
        std_ok=nanstd(x0s_ok(:,3));
        var_ok=nanvar(x0s_ok(:,3));
        CV_ok=std_ok/mean_ok;
        MSPE_ok=sqrt(nanmean((R(:,3)-x0s_ok(:,3)).^2));
        CC_ok=corrcoef(R(:,3),x0s_ok(:,3),'rows','pairwise');
        
        %CK
        itype2=1.7;
        [x0s_ck,~,~,~,~,~,~,nbimag_ck]=cokrioct(x,G,model,c,itype2,avg,block,nd,ival,nk_ck,rad,ntok);
        nbimag_ck=sum(nbimag_ck);
        neg_ck=length(find(x0s_ck(:,3)<0));
        nan_ck=length(find(isnan(x0s_ck)));
        
        x0s_ck(x0s_ck<0)=0;
        mean_ck=nanmean(x0s_ck(:,3));
        std_ck=nanstd(x0s_ck(:,3));
        var_ck=nanvar(x0s_ck(:,3));
        CV_ck=std_ck/mean_ck;
        MSPE_ck=sqrt(nanmean((R(:,3)-x0s_ck(:,3)).^2));
        CC_ck=corrcoef(R(:,3),x0s_ck(:,3),'rows','pairwise');
        
        %UC
        y=x(:,3);
        [z]=anamor(y);
        nscore_WL=[z(:,2),z(:,1)];
        
        phi=nscore_WL;
        phi_s=sortrows(phi,1);
        phi_s=[-7 0;phi_s;7 1.001*max(phi(:,2))];
        u=norminv([0.0005:0.001:1]',0,1);
        z=interp1(phi_s(:,1),phi_s(:,2),u,'linear',phi_s(end,2));
        z(1)=0;
        phi=[-7 0;u z;7 1.001*z(end)];
        itype3=3;
        [x0s_uc,~,~,~,~,~,~,~]=cokrioct(x,G_UC,model,c,itype3,avg,panneau,nd_uc,ival,nk_uc,rad,ntok);
        
        zk=x0s_uc(:,3);
        
        [T_uc,m_uc,~,~,~,~,~,~,~]=cond_unif_v3(zk,model,c,phi,panneau,smu,vc);
        
        ton_uc=mean(T_uc);
        mean_uc=mean(m_uc);
        
        
    elseif dec==1
        %OK
        itype=3;
        [x0s_ok,~,~,~,~,~,~,nbimag_ok]=cokrioct(x,G,model,c,itype,avg,block,nd,ival,nk_ok,rad,ntok);
        nbimag_ok=sum(nbimag_ok);
        neg_ok=length(find(x0s_ok(:,3)<0));
        nan_ok=length(find(isnan(x0s_ok)));
        
        x0s_ok(x0s_ok<0)=0;
        mean_ok=nanmean(x0s_ok(:,3));
        std_ok=nanstd(x0s_ok(:,3));
        var_ok=nanvar(x0s_ok(:,3));
        CV_ok=std_ok/mean_ok;
        MSPE_ok=sqrt(nanmean((R(:,3)-x0s_ok(:,3)).^2));
        CC_ok=corrcoef(R(:,3),x0s_ok(:,3),'rows','pairwise');
        
        %CK
        itype2=1.7;
        [x0s_ck,~,~,~,~,~,~,nbimag_ck]=cokrioct(x,G,model,c,itype2,avg,block,nd,ival,nk_ck,rad,ntok);
        nbimag_ck=sum(nbimag_ck);
        neg_ck=length(find(x0s_ck(:,3)<0));
        nan_ck=length(find(isnan(x0s_ck)));
        
        x0s_ck(x0s_ck<0)=0;
        mean_ck=nanmean(x0s_ck(:,3));
        std_ck=nanstd(x0s_ck(:,3));
        var_ck=nanvar(x0s_ck(:,3));
        CV_ck=std_ck/mean_ck;
        MSPE_ck=sqrt(nanmean((R(:,3)-x0s_ck(:,3)).^2));
        CC_ck=corrcoef(R(:,3),x0s_ck(:,3),'rows','pairwise');
        
        %UC
        y=x_dec(:,3);
        [z]=anamor(y);
        nscore_WL=[z(:,2),z(:,1)];
        
        phi=nscore_WL;
        phi_s=sortrows(phi,1);
        phi_s=[-7 0;phi_s;7 1.001*max(phi(:,2))];
        u=norminv([0.0005:0.001:1]',0,1);
        z=interp1(phi_s(:,1),phi_s(:,2),u,'linear',phi_s(end,2));
        %z(1)=0;
        phi=[-7 0;u z;7 1.001*z(end)];
        itype3=3;
        [x0s_uc,~,~,~,~,~,~,~]=cokrioct(x,G_UC,model,c,itype3,avg,panneau,nd_uc,ival,nk_uc,rad,ntok);
        
        zk=x0s_uc(:,3);
        
        [T_uc,m_uc,~,~,~,~,~,~,~]=cond_unif_v3(zk,model,c,phi,panneau,smu,vc);
        
        ton_uc=mean(T_uc);
        mean_uc=mean(m_uc);
    end
    % Recovery functions
    OK=x0s_ok(:,3);
    CK=x0s_ck(:,3);
    RR=R(:,3);
    c=vc;
    
    for i=1:length(c)
        tc_r(i)=mean(RR>=c(i));
        tc_ok(i)=mean(OK>=c(i));
        tc_ck(i)=mean(CK>=c(i));
        
        mc_r(i)=mean(RR(RR>=c(i)));
        mc_ok(i)=mean(OK(OK>=c(i)));
        mc_ck(i)=mean(CK(CK>=c(i)));
    end
    
    % Stats
    stat.nan_ok=nan_ok;
    stat.nan_ck=nan_ck;
    stat.nbneg_ok=neg_ok;
    stat.nbneg_ck=neg_ck;
    stat.nbimag_ok=nbimag_ok;
    stat.nbimag_ck=nbimag_ck;
    stat.osr_ok=100*((var(R(:,3))-nanvar(x0s_ok(:,3)))/var(R(:,3)));
    stat.osr_ck=100*((var(R(:,3))-nanvar(x0s_ck(:,3)))/var(R(:,3)));
    stat.mspe_ok=MSPE_ok;
    stat.mspe_ck=MSPE_ck;
    stat.corr_ok=CC_ok(1,2);
    stat.corr_ck=CC_ck(1,2);
    stat.tc_r=tc_r;
    stat.mc_r=mc_r;
    stat.tc_ok=tc_ok;
    stat.mc_ok=mc_ok;
    stat.tc_ck=tc_ck;
    stat.mc_ck=mc_ck;
    stat.summary=[nan_ck,neg_ck,nbimag_ck];
    stat.OK=[length(x0s_ok(:,3)),mean_ok,var_ok,CV_ok,MSPE_ok,CC_ok(1,2)];
    stat.CK=[length(x0s_ck(:,3)),mean_ck,var_ck,CV_ck,MSPE_ck,CC_ck(1,2)];
    
    % Plot
    figure(12+cas);clf
    plot(c,tc_r,'--r',c,tc_ok,'-b',c,tc_ck,'-k','linewidth',2)
    hold on
    plot(c,ton_uc,'-C','linewidth',2)
    grid on
    h=legend('Real','OK','CK','UC');
    set(h,'Interpreter','Latex','Location','northeast','FontSize',12);
    set(gca,'TickLabelInterpreter','latex','FontSize',12)
    box off
    xlabel('Cut-off','Interpreter','Latex','FontSize', 12)
    ylabel('Tonnage','Interpreter','Latex','FontSize', 12)
    %title(['Curve tonnage - cut-off: Neighbors ',num2str(m)],'Interpreter','Latex','FontSize', 12)
    print(12+cas,'-depsc2',['Ton_WL_cas' num2str(cas) '.eps']);
    close(figure(12+cas));
    
    
    figure(13+cas);clf
    plot(c,mc_r,'--r',c,mc_ok,'-b',c,mc_ck,'-k','linewidth',2)
    hold on
    plot(c,mean_uc,'-c','linewidth',2)
    grid on
    h=legend('Real','OK','CK','UC');
    set(h,'Interpreter','Latex','Location','southeast','FontSize',12);
    set(gca,'TickLabelInterpreter','latex','FontSize',12)
    box off
    xlabel('Cut-off','Interpreter','Latex','FontSize', 12)
    ylabel('Ore grade','Interpreter','Latex','FontSize', 12)
    %title(['Curve tonnage - cut-off: Neighbors ',num2str(m)],'Interpreter','Latex','FontSize', 12)
    print(13+cas,'-depsc2',['Ore_WL_cas' num2str(cas) '.eps']);
    close(figure(13+cas));
    
elseif x2==4
    
    if dec==0
        %    OK
        itype=3;
        [x0s_ok,~,~,~,~,~,~,nbimag_ok]=cokrioct(x,G,model,c,itype,avg,block,nd,ival,nk_ok,rad,ntok);
        nbimag_ok=sum(nbimag_ok);
        neg_ok=length(find(x0s_ok(:,4)<0));
        nan_ok=length(find(isnan(x0s_ok)));
        
        x0s_ok(x0s_ok<0)=0;
        mean_ok=nanmean(x0s_ok(:,4));
        std_ok=nanstd(x0s_ok(:,4));
        var_ok=nanvar(x0s_ok(:,4));
        CV_ok=std_ok/mean_ok;
        MSPE_ok=sqrt(nanmean((R(:,4)-x0s_ok(:,4)).^2));
        CC_ok=corrcoef(R(:,4),x0s_ok(:,4),'rows','pairwise');
        
        %CK
        itype2=1.7;
        [x0s_ck,~,~,~,~,~,~,nbimag_ck]=cokrioct(x,G,model,c,itype2,avg,block,nd,ival,nk_ck,rad,ntok);
        nbimag_ck=sum(nbimag_ck);
        neg_ck=length(find(x0s_ck(:,4)<0));
        nan_ck=length(find(isnan(x0s_ck)));
        
        x0s_ck(x0s_ck<0)=0;
        mean_ck=nanmean(x0s_ck(:,4));
        std_ck=nanstd(x0s_ck(:,4));
        var_ck=nanvar(x0s_ck(:,4));
        CV_ck=std_ck/mean_ck;
        MSPE_ck=sqrt(nanmean((R(:,4)-x0s_ck(:,4)).^2));
        CC_ck=corrcoef(R(:,4),x0s_ck(:,4),'rows','pairwise');
        
        %UC
        
        itype3=3;
        [x0s_uc,~,~,~,~,~,~,~]=cokrioct(x,G_UC,model,c,itype3,avg,panneau,nd_uc,ival,nk_uc,rad,ntok);
        
        y=x(:,4);
        [z]=anamor(y);
        nscore_Gold=[z(:,2),z(:,1)];
        
        phi=nscore_Gold;
        phi_s=sortrows(phi,1);
        phi_s=[-7 0;phi_s;7 1.001*max(phi(:,2))];
        u=norminv([0.0000000001:0.001:1]',0,1);
        z=interp1(phi_s(:,1),phi_s(:,2),u,'linear',phi_s(end,2));
        z(1)=0;
        phi=[-7 0;u z;7 1.001*z(end)];
        
        zk=x0s_uc(:,4);
        
        [T_uc,m_uc,~,~,~,~,~,~,~]=cond_unif_v3(zk,model,c,phi,panneau,smu,vc);
        
        ton_uc=mean(T_uc);
        mean_uc=mean(m_uc);
        
    elseif dec==1
        
        %OK
        %    OK
        itype=3;
        [x0s_ok,~,~,~,~,~,~,nbimag_ok]=cokrioct(x,G,model,c,itype,avg,block,nd,ival,nk_ok,rad,ntok);
        nbimag_ok=sum(nbimag_ok);
        neg_ok=length(find(x0s_ok(:,4)<0));
        nan_ok=length(find(isnan(x0s_ok)));
        
        x0s_ok(x0s_ok<0)=0;
        mean_ok=nanmean(x0s_ok(:,4));
        std_ok=nanstd(x0s_ok(:,4));
        var_ok=nanvar(x0s_ok(:,4));
        CV_ok=std_ok/mean_ok;
        MSPE_ok=sqrt(nanmean((R(:,4)-x0s_ok(:,4)).^2));
        CC_ok=corrcoef(R(:,4),x0s_ok(:,4),'rows','pairwise');
        
        %CK
        itype2=1.7;
        [x0s_ck,~,~,~,~,~,~,nbimag_ck]=cokrioct(x,G,model,c,itype2,avg,block,nd,ival,nk_ck,rad,ntok);
        nbimag_ck=sum(nbimag_ck);
        neg_ck=length(find(x0s_ck(:,4)<0));
        nan_ck=length(find(isnan(x0s_ck)));
        
        x0s_ck(x0s_ck<0)=0;
        mean_ck=nanmean(x0s_ck(:,4));
        std_ck=nanstd(x0s_ck(:,4));
        var_ck=nanvar(x0s_ck(:,4));
        CV_ck=std_ck/mean_ck;
        MSPE_ck=sqrt(nanmean((R(:,4)-x0s_ck(:,4)).^2));
        CC_ck=corrcoef(R(:,4),x0s_ck(:,4),'rows','pairwise');
        
        %UC
        y=x_dec(:,4);
        [z]=anamor(y);
        nscore_Gold=[z(:,2),z(:,1)];
        
        phi=nscore_Gold;
        phi_s=sortrows(phi,1);
        phi_s=[-7 0;phi_s;7 1.001*max(phi(:,2))];
        u=norminv([0.0005:0.001:1]',0,1);
        z=interp1(phi_s(:,1),phi_s(:,2),u,'linear',phi_s(end,2));
        z(1)=0;
        phi=[-7 0;u z;7 1.001*z(end)];
        itype3=3;
        [x0s_uc,~,~,~,~,~,~,~,~]=cokri(x,G_UC,model,c,itype3,avg,panneau,nd_uc,ival,nk_uc,rad,ntok);
        
        zk=x0s_uc(:,4);
        
        [T_uc,m_uc,~,~,~,~,~,~,~]=cond_unif_v3(zk,model,c,phi,panneau,smu,vc);
        
        ton_uc=mean(T_uc);
        mean_uc=mean(m_uc);
        
    end
    
    % Recovery functions
    OK=x0s_ok(:,4);
    CK=x0s_ck(:,4);
    RR=R(:,4);
    c=vc;
    
    for i=1:length(c)
        tc_r(i)=mean(RR>=c(i));
        tc_ok(i)=mean(OK>=c(i));
        tc_ck(i)=mean(CK>=c(i));
        
        mc_r(i)=mean(RR(RR>=c(i)));
        mc_ok(i)=mean(OK(OK>=c(i)));
        mc_ck(i)=mean(CK(CK>=c(i)));
    end
    
    % Stats
    stat.nan_ok=nan_ok;
    stat.nan_ck=nan_ck;
    stat.nbneg_ok=neg_ok;
    stat.nbneg_ck=neg_ck;
    stat.nbimag_ok=nbimag_ok;
    stat.nbimag_ck=nbimag_ck;
    stat.osr_ok=100*((var(R(:,4))-var(x0s_ok(:,4)))/var(R(:,4)));
    stat.osr_ck=100*((nanvar(R(:,4))-nanvar(x0s_ck(:,4)))/nanvar(R(:,4)));
    stat.mspe_ok=MSPE_ok;
    stat.mspe_ck=MSPE_ck;
    stat.corr_ok=CC_ok(1,2);
    stat.corr_ck=CC_ck(1,2);
    stat.tc_r=tc_r;
    stat.mc_r=mc_r;
    stat.tc_ok=tc_ok;
    stat.mc_ok=mc_ok;
    stat.tc_ck=tc_ck;
    stat.mc_ck=mc_ck;
    stat.summary=[nan_ck,neg_ck,nbimag_ck];
    stat.OK=[length(x0s_ok(:,4)),mean_ok,var_ok,CV_ok,MSPE_ok,CC_ok(1,2)];
    stat.CK=[length(x0s_ck(:,4)),mean_ck,var_ck,CV_ck,MSPE_ck,CC_ck(1,2)];
    
    %        Plot
    figure(12+cas);clf
    plot(c,tc_r,'--r',c,tc_ok,'-b',c,tc_ck,'-k','linewidth',2)
    hold on
    plot(c,ton_uc,'-C','linewidth',2)
    grid on
    h=legend('Real','OK','CK','UC');
    set(h,'Interpreter','Latex','Location','northeast','FontSize',12);
    set(gca,'TickLabelInterpreter','latex','FontSize',12)
    box off
    xlabel('Cut-off','Interpreter','Latex','FontSize', 12)
    ylabel('Tonnage','Interpreter','Latex','FontSize', 12)
    %title(['Curve tonnage - cut-off: Neighbors ',num2str(m)],'Interpreter','Latex','FontSize', 12)
    print(12+cas,'-depsc2',['Ton_Gold_cas' num2str(cas) '.eps']);
    close(figure(12+cas));
    
    figure(13+cas);clf
    plot(c,mc_r,'--r',c,mc_ok,'-b',c,mc_ck,'-k','linewidth',2)
    hold on
    plot(c,mean_uc,'-c','linewidth',2)
    grid on
    h=legend('Real','OK','CK','UC');
    set(h,'Interpreter','Latex','Location','southeast','FontSize',12);
    set(gca,'TickLabelInterpreter','latex','FontSize',12)
    box off
    xlabel('Cut-off','Interpreter','Latex','FontSize', 12)
    ylabel('Ore grade','Interpreter','Latex','FontSize', 12)
    %title(['Curve tonnage - cut-off: Neighbors ',num2str(m)],'Interpreter','Latex','FontSize', 12)
    print(13+cas,'-depsc2',['Ore_Gold_cas' num2str(cas) '.eps']);
    close(figure(13+cas));
    drawnow
end

end
