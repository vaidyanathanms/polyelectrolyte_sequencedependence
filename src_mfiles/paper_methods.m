%% To compute the adsorbed values in different formats

clear
close all
clc
format long

%% Plot Colors

green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17];
pclr = {'r','b',green,'k','m', gold, orange};
lsty = {'-','--',':'};
msty = {'d','s','o','x'};

%% Flags

fig1a = 0;
fig1b = 1;
fig2a = 1;
fig3a = 0;
figr1 = 0;

%% Inputs

nmonfree = 30; nmongraft = 30; ngraft = 64;nbase = 10*ngraft;nsalt = 510;
nfreearr = [100];
ncntr = abs((nfreearr*nmonfree-ngraft*nmongraft)/2);
ntotarr = nfreearr*nmonfree + nmongraft*ngraft + nbase + nsalt + ncntr;
cutoff = 0.9; lz = 120; area=53^2;
rhofree = nfreearr*30/(lz*area);
rcut = '1.10';
nbins = 120;

%% Figure 1(a) - Ratio of Free-Anions within cutoff any brush

if fig1a ~= 0
    
    for nvals = 1:length(nfreearr)
        
        nfree = nfreearr(nvals);
        fnr = fopen(sprintf('./All_TxtFiles_Paper/fig1a_BoundRat_%d_%s.txt',nfreearr(nvals),rcut),'w');
        fprintf(fnr,'%s \n','Method: Explicitly Compute anions within,rcut=1.1sigma)');
        fprintf(fnr,'%s\t %s\n','Arch','f(r_cut)');
        
        for i = 1:4
            
            if i == 1
                dirstr = 'bl_bl';
            elseif i == 2
                dirstr = 'bl_al';
            elseif i == 3
                dirstr = 'al_bl';
            else
                dirstr = 'al_al';
            end
            
            fid = fopen(sprintf('./n_%d/results_%d_%s/adsfracv2_rcut_%s.lammpstrj',nfree,nfree,dirstr,rcut));
            data = textscan(fid,'%f%f%f');
            fld  = cell2mat(data);
            frac = fld(:,2);
            fracmean = mean(frac);
            fprintf(fnr,'%s\t%g\n',dirstr,fracmean);
            fclose(fid);
            
        end
        
        fclose(fnr);
        
    end
    
    
    % Plot the frac values - Fig 1(a) Ratio of Anion to Cation
    
    data_bb = zeros(length(nfreearr),1);
    data_ab = zeros(length(nfreearr),1);
    data_ba = zeros(length(nfreearr),1);
    data_aa = zeros(length(nfreearr),1);
    
    for i = 1:length(nfreearr)
        
        fnr = fopen(sprintf('./All_TxtFiles_Paper/fig1a_BoundRat_%d_%s.txt',nfreearr(i),rcut),'r');
        data = textscan(fnr,'%s%f','Headerlines',2);
        
        fld = cell2mat(data(2));
        data_bb(i,1) = fld(1,1);
        data_ba(i,1) = fld(2,1);
        data_ab(i,1) = fld(3,1);
        data_aa(i,1) = fld(4,1);
        
        fclose(fnr);
        
    end
    
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$N_p/N_g$','FontSize',20,'Interpreter','Latex')
    ylabel('$f$','FontSize',20,'Interpreter','Latex')
    
    plot(nfreearr/ngraft,data_bb/(ngraft*nmongraft),'color',pclr{1},'LineWidth',2,'LineStyle',lsty{3},'Marker',msty{1},'MarkerSize',8,'MarkerFaceColor',pclr{1})
    plot(nfreearr/ngraft,data_ba/(ngraft*nmongraft),'color',pclr{2},'LineWidth',2,'LineStyle',lsty{3},'Marker',msty{2},'MarkerSize',8,'MarkerFaceColor',pclr{2})
    plot(nfreearr/ngraft,data_ab/(ngraft*nmongraft),'color',pclr{3},'LineWidth',2,'LineStyle',lsty{3},'Marker',msty{3},'MarkerSize',8,'MarkerFaceColor',pclr{3})
    plot(nfreearr/ngraft,data_aa/(ngraft*nmongraft),'color',pclr{4},'LineWidth',2,'LineStyle',lsty{3},'Marker',msty{4},'MarkerSize',8,'MarkerFaceColor',pclr{4})
    
    
    legendinfo{1} = 'Block-Block';
    legendinfo{2} = 'Block-Alter';
    legendinfo{3} = 'Alter-Block';
    legendinfo{4} = 'Alter-Alter';
    
    legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
    legend boxoff
    saveas(h1,sprintf('Fig1a_CountAnions_%s',rcut),'png');
    clear legendinfo
    
end

%% Figure 1(b) - Net Charge withing the brush regime

if fig1b ~= 0
    
    chargearr = [0;1;0;-1;1;-1];
    
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$N_{pa}/N_{pc}$','FontSize',20,'Interpreter','Latex')
    ylabel('$f$','FontSize',20,'Interpreter','Latex')
    
    for nvals = 1:length(nfreearr)
        
        nfree = nfreearr(nvals);
        fnr = fopen(sprintf('./fig1b_QnetBound_%d.txt',nfreearr(nvals)),'w');
        fprintf(fnr,'%s \n','NetCharge: Q_{b}=\Delta(n_g)*\int(\sum(q_j n_j(z)dz, j=all entities)z=0,Lz))');
        fprintf(fnr,'%s\t %s\n','Arch','Q_{b}');
        
        if nfree > ngraft
            cntrcharg = 1;
        elseif nfree == ngraft
            cntrcharg = 0;
        elseif nfree < ngraft
            cntrcharg = -1;
        end
        
        for i = 1:4
            
            if i == 1
                dirstr = 'bl_bl';
            elseif i == 2
                dirstr = 'bl_al';
            elseif i == 3
                dirstr = 'al_bl';
            elseif i == 4
                dirstr = 'al_al';
            else
                disp('No Correct String')
                break;
            end
            
            
            fid = fopen(sprintf('./results_dens/results_%d_%s/dens.lammpstrj',nfree,dirstr));
            data = textscan(fid,'%f%f%f%f%f%f%f%f','Headerlines',1);
            fld = cell2mat(data);
            rdata     = fld(:,1);
            pneutral_brush = fld(:,2)*ngraft*nmongraft/2*lz*area/nbins;
            pcation_brush = fld(:,3)*ngraft*nmongraft/2*lz*area/nbins;
            pneutral_free  = fld(:,4)*nfreearr(nvals)*nmonfree/2*lz*area/nbins;
            panion_brush = fld(:,5)*nfreearr(nvals)*nmonfree/2*lz*area/nbins;
            pposions  = fld(:,6)*nsalt*lz*area/nbins;
            pnegions  = fld(:,7)*nsalt*lz*area/nbins;
            pcntrions = fld(:,8)*ncntr(nvals)*lz*area/nbins;
            
            if nfree*nmonfree ~= ngraft*nmongraft
                netcharge = chargearr(1)*pneutral_brush + chargearr(2)*pcation_brush + ...
                    chargearr(3)*pneutral_free + chargearr(4)*panion_brush + ...
                    chargearr(5)*pposions + chargearr(6)*pnegions + cntrcharg*pcntrions;
            else
                netcharge = chargearr(1)*pneutral_brush + chargearr(2)*pcation_brush + ...
                    chargearr(3)*pneutral_free + chargearr(4)*panion_brush + ...
                    chargearr(5)*pposions + chargearr(6)*pnegions;
            end
            
            sumq = 0.5*netcharge(1);
            
            fid_g  = fopen(sprintf('./results_dens/results_%d_%s/grpdens.lammpstrj',nfree,dirstr));
            data_g = textscan(fid_g,'%f%f%f','Headerlines',1);
            
            fld_g   = cell2mat(data_g);
            dens_g  = fld_g(:,2);
            [maxden, imaxden]  = max(dens_g);
            
            trapz(rdata,pcntrions)
            
            % Find edge of brush
            
            for k = imaxden:length(dens_g)
                
                if dens_g(k,1) < 0.05*maxden
                    
                    i_edge = k;
                    break;
                    
                end
                
            end
            
            qofr = zeros(i_edge,1); rqofr = zeros(i_edge,1);
            qofr(1,1) = sumq; rqofr(1,1) = 0.5*rdata(1);
            
            % Integrate charge to the edge of the brush
            for k = 2:i_edge
                
                sumq = sumq + 0.5*(rdata(k)-rdata(k-1))*(netcharge(k)+netcharge(k-1));
                qofr(k,1)  = sumq;
                rqofr(k,1) = 0.5*(rdata(k)+rdata(k+1));
                
            end
            
            fprintf(fnr,'%s\t%g\n',dirstr,sumq);
            fclose(fid_g);
            fclose(fid);
            
            if(nvals == 1)
                
                plot(rqofr,qofr,'color',pclr{i},'LineWidth',2,'LineStyle',lsty{1})
                
            end
            
        end
        
        fclose(fnr);
        
    end
    
    % Plot Q_{b} data
    
    data_bb = zeros(length(nfreearr),1);
    data_ab = zeros(length(nfreearr),1);
    data_ba = zeros(length(nfreearr),1);
    data_aa = zeros(length(nfreearr),1);
    
    for i = 1:length(nfreearr)
        
        fnr = fopen(sprintf('./fig1b_QnetBound_%d.txt',nfreearr(i)),'r');
        data = textscan(fnr,'%s%f','Headerlines',2);
        
        fld = cell2mat(data(2));
        data_bb(i,1) = fld(1,1);
        data_ba(i,1) = fld(2,1);
        data_ab(i,1) = fld(3,1);
        data_aa(i,1) = fld(4,1);
        
        fclose(fnr);
        
    end
    
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$N_{pa}/N_{pc}$','FontSize',20,'Interpreter','Latex')
    ylabel('$Q_{b}$','FontSize',20,'Interpreter','Latex')
    
    plot(nfreearr/ngraft,data_bb,'color',pclr{1},'LineWidth',2,'LineStyle',lsty{3},'Marker',msty{1},'MarkerSize',8,'MarkerFaceColor',pclr{1})
    plot(nfreearr/ngraft,data_ba,'color',pclr{2},'LineWidth',2,'LineStyle',lsty{3},'Marker',msty{2},'MarkerSize',8,'MarkerFaceColor',pclr{2})
    plot(nfreearr/ngraft,data_ab,'color',pclr{3},'LineWidth',2,'LineStyle',lsty{3},'Marker',msty{3},'MarkerSize',8,'MarkerFaceColor',pclr{3})
    plot(nfreearr/ngraft,data_aa,'color',pclr{4},'LineWidth',2,'LineStyle',lsty{3},'Marker',msty{4},'MarkerSize',8,'MarkerFaceColor',pclr{4})
    
    
    legendinfo{1} = 'Block-Block';
    legendinfo{2} = 'Block-Alter';
    legendinfo{3} = 'Alter-Block';
    legendinfo{4} = 'Alter-Alter';
    
    legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
    legend boxoff
    saveas(h1,'Fig1b_QnetBound_%s','png');
    clear legendinfo
    
end

%% Figure 2a - q(z) vs z

if fig2a ~= 0
    
    chargearr = [0;1;0;-1;1;-1];
    
    for nvals = 1:length(nfreearr)
        
        nfree = nfreearr(nvals);
        
        
        if nfree > ngraft
            cntrcharg = 1;
        elseif nfree == ngraft
            cntrcharg = 0;
        elseif nfree < ngraft
            cntrcharg = -1;
        end
        
        for i = 1:4
            
            if i == 1
                dirstr = 'bl_bl';
            elseif i == 2
                dirstr = 'bl_al';
            elseif i == 3
                dirstr = 'al_bl';
            elseif i == 4
                dirstr = 'al_al';
            else
                disp('No Correct String')
                break;
            end
            
            fnr = fopen(sprintf('./fig2a_Qz_%d_%s.txt',nfreearr(nvals),dirstr),'w');
            fprintf(fnr,'%s \n','NetCharge: Q_{z}=\int(\sum(q_j n_j(z)dz, j=all entities)z=0,z))');
            fprintf(fnr,'%s\t %s\n','z','q_{z}');
            
            fid = fopen(sprintf('./results_dens/results_%d_%s/dens.lammpstrj',nfree,dirstr));
            data = textscan(fid,'%f%f%f%f%f%f%f%f','Headerlines',1);
            fld = cell2mat(data);
            rdata     = fld(:,1);
            pneutral_brush = fld(:,2)*ngraft*nmongraft/2*lz*area/nbins;
            pcation_brush = fld(:,3)*ngraft*nmongraft/2*lz*area/nbins;
            pneutral_free  = fld(:,4)*nfreearr(nvals)*nmonfree/2*lz*area/nbins;
            panion_brush = fld(:,5)*nfreearr(nvals)*nmonfree/2*lz*area/nbins;
            pposions  = fld(:,6)*nsalt*lz*area/nbins;
            pnegions  = fld(:,7)*nsalt*lz*area/nbins;
            pcntrions = fld(:,8)*ncntr(nvals)*lz*area/nbins;
            
            if nfree*nmonfree ~= ngraft*nmongraft
                netcharge = chargearr(1)*pneutral_brush + chargearr(2)*pcation_brush + ...
                    chargearr(3)*pneutral_free + chargearr(4)*panion_brush + ...
                    chargearr(5)*pposions + chargearr(6)*pnegions + cntrcharg*pcntrions;
            else
                netcharge = chargearr(1)*pneutral_brush + chargearr(2)*pcation_brush + ...
                    chargearr(3)*pneutral_free + chargearr(4)*panion_brush + ...
                    chargearr(5)*pposions + chargearr(6)*pnegions;
            end
            
            sumq = 0.5*netcharge(1);
            
            fid_g  = fopen(sprintf('./results_dens/results_%d_%s/grpdens.lammpstrj',nfree,dirstr));
            data_g = textscan(fid_g,'%f%f%f','Headerlines',1);
            
            fld_g   = cell2mat(data_g);
            dens_g  = fld_g(:,2);
            [maxden, imaxden]  = max(dens_g);
            
            qofr = zeros(length(netcharge),1); rqofr = zeros(length(netcharge),1);
            qofr(1,1) = sumq; rqofr(1,1) = 0.5*rdata(1);
            
            % Integrate charge to the edge of the brush
            for k = 2:length(netcharge)
                
                sumq = sumq + 0.5*(rdata(k)-rdata(k-1))*(netcharge(k)+netcharge(k-1));
                qofr(k,1)  = sumq;
                rqofr(k,1) = 0.5*(rdata(k)+rdata(k-1));
                
            end
            
            for k = 1:length(netcharge)
                
                fprintf(fnr,'%g\t%g\n',rqofr(k,1),qofr(k,1));
                
            end
            
            fclose(fid_g);
            fclose(fid);
            fclose(fnr);
            
        end
        
    end
    
    % Plot Q_{z} data
    
    for i = 1:length(nfreearr)
        
        h1 = figure;
        hold on
        box on
        set(gca,'FontSize',16)
        xlabel('$z/L_z$','FontSize',20,'Interpreter','Latex')
        ylabel('$q(z)$','FontSize',20,'Interpreter','Latex')
        
        for j = 1:4
            
            if j == 1
                dirstr = 'bl_bl';
            elseif j == 2
                dirstr = 'bl_al';
            elseif j == 3
                dirstr = 'al_bl';
            elseif j == 4
                dirstr = 'al_al';
            else
                disp('No Correct String')
                break;
            end
            
            fnr = fopen(sprintf('./fig2a_Qz_%d_%s.txt',nfreearr(i),dirstr),'r');
            data = textscan(fnr,'%f%f','Headerlines',2);
            
            pl_data = zeros(length(rqofr),2);
            fld = cell2mat(data);
            pl_data(:,1) = fld(:,1);
            pl_data(:,2) = fld(:,2);
            plot(pl_data(:,1)/lz,pl_data(:,2),'color',pclr{j},'LineWidth',2,'LineStyle',lsty{1})
            fclose(fnr);
            
        end
        
        legendinfo{1} = 'Block-Block';
        legendinfo{2} = 'Block-Alter';
        legendinfo{3} = 'Alter-Block';
        legendinfo{4} = 'Alter-Alter';
        
        legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
        legend boxoff
        saveas(h1,sprintf('Fig2a_Qz_%g',nfreearr(i)),'png');
        clear legendinfo
        
    end
    
end

%% Figure 3a - density profiles

if fig3a ~= 0
    
    for nvals = 1:length(nfreearr)
        
        nfree = nfreearr(nvals);
        ntotmons = ntotarr(nvals);
        
        h1 = figure;
        hold on
        box on
        set(gca,'FontSize',16)
        h = zeros(4,1);
        for i = 1:4
            
            h(i) = subplot(4,1,i);
            
            if i == 1
                dirstr = 'bl_bl';
            elseif i == 2
                dirstr = 'bl_al';
            elseif i == 3
                dirstr = 'al_bl';
            elseif i == 4
                dirstr = 'al_al';
            else
                disp('No Correct String')
                break;
            end
            
            fid = fopen(sprintf('./n_%d/results_%d_%s/grpdens.lammpstrj',nfree,nfree,dirstr));
            data = textscan(fid,'%f%f%f','Headerlines',1);
            fld = cell2mat(data);
            
            rdata     = fld(:,1);
            pl_data(:,1) = fld(:,2)*ngraft*nmongraft/ntotmons;
            pl_data(:,2) = fld(:,3)*nfree*nmonfree/ntotmons;
            
            patch(rdata/lz,pl_data(:,1),'r');
            alpha(0.3);
            patch(rdata/lz,pl_data(:,2),'g')
            alpha(0.5);
            xlim([0 0.3])
            fclose(fid);
            
        end
        
        set(h(1),'xticklabel',[]);
        set(h(2),'xticklabel',[]);
        set(h(3),'xticklabel',[]);
        
        pos=get(h,'position');
        bottom=pos{4}(2);
        top=pos{1}(2)+pos{1}(4);
        plotspace=top-bottom;
        
        pos{4}(4)=plotspace/4;
        pos{3}(4)=plotspace/4;
        pos{2}(4)=plotspace/4;
        pos{1}(4)=plotspace/4;
        
        pos{1}(2)=bottom + 3*plotspace/4;
        pos{2}(2)=bottom + 2*plotspace/4;
        pos{3}(2)=bottom + plotspace/4;
        
        
        set(h(1),'position',pos{1},'FontSize',16);
        set(h(2),'position',pos{2},'FontSize',16);
        set(h(3),'position',pos{3},'FontSize',16);
        set(h(4),'position',pos{4},'FontSize',16);
        box on
        
        yticks(h(4),[])
        if(nfree >= 80)
            yticks(h(1),[0 2.5e-6 5e-6 7.5e-6])
            yticklabels(h(1),[0 0.25 0.5 0.75])
            yticks(h(2),[0 2.5e-6 5e-6 7.5e-6])
            yticklabels(h(2),[0 0.25 0.5 0.75])
            yticks(h(3),[0 2.5e-6 5e-6 7.5e-6])
            yticklabels(h(3),[0 0.25 0.5 0.75])
            yticks(h(4),[0 2.5e-6 5e-6 7.5e-6])
            yticklabels(h(4),[0 0.25 0.5 0.75])
        else
            yticks(h(1),[0 5e-6 1e-5 1.5e-5])
            yticklabels(h(1),[0 0.5 1.0 1.5])
            yticks(h(2),[0 5e-6 1e-5 1.5e-5])
            yticklabels(h(2),[0 0.5 1.0 1.5])
            yticks(h(3),[0 5e-6 1e-5 1.5e-5])
            yticklabels(h(3),[0 0.5 1.0 1.5])
            yticks(h(4),[0 5e-6 1e-5 1.5e-5])
            yticklabels(h(4),[0 0.5 1.0 1.5])
            
        end
            
        xlabel('$z/L_z$','FontSize',20,'Interpreter','Latex')
        ylabel('$\rho(z) \times 10^5$','FontSize',20,'Interpreter','Latex')
        box(h(1),'on');
        box(h(2),'on');
        box(h(3),'on');
        legendinfo{1} = '$\rho_g(z)$';
        legendinfo{2} = '$\rho_p(z)$';
        legend(h(1), legendinfo, 'Interpreter','Latex','FontSize',16)
        box(legend,'off')
    end
    
end

%% Rest of the Figures - 1

if figr1 ~= 0
    
    chargearr = [0;1;0;-1;1;-1];
    nfreecnt  = 0;
    for nvals = 1:length(nfreearr)
       
        
        nfree = nfreearr(nvals);
        
%         if nfree >= ngraft
%             
%             break
%             
%         else
            
            nfreecnt = nfreecnt + 1;
            
%         end
        
        if nfree > ngraft
            cntrcharg = 1;
        elseif nfree == ngraft
            cntrcharg = 0;
        elseif nfree < ngraft
            cntrcharg = -1;
        end
        
        for i = 1:4
            
            if i == 1
                dirstr = 'bl_bl';
            elseif i == 2
                dirstr = 'bl_al';
            elseif i == 3
                dirstr = 'al_bl';
            elseif i == 4
                dirstr = 'al_al';
            else
                disp('No Correct String')
                break;
            end
            
            fnr = fopen(sprintf('./All_TxtFiles_Paper/figr1_ionQz_%d_%s.txt',nfreearr(nvals),dirstr),'w');
            fprintf(fnr,'%s \n','NetCharge: Q_{z}=\int(\sum(q_j n_j(z)dz, j=all ions)z=0,z))');
            fprintf(fnr,'%s\t %s\t %s\n','z','qpos_{z}','qneg_{z}');
            
            fid = fopen(sprintf('./n_%d/results_%d_%s/dens.lammpstrj',nfree,nfree,dirstr));
            data = textscan(fid,'%f%f%f%f%f%f%f%f','Headerlines',1);
            fld = cell2mat(data);
            rdata     = fld(:,1);
            pneutral_brush = fld(:,2)*ngraft*nmongraft/2*lz*area/nbins;
            pcation_brush = fld(:,3)*ngraft*nmongraft/2*lz*area/nbins;
            pneutral_free  = fld(:,4)*nfreearr(nvals)*nmonfree/2*lz*area/nbins;
            panion_brush = fld(:,5)*nfreearr(nvals)*nmonfree/2*lz*area/nbins;
            pposions  = fld(:,6)*nsalt*lz*area/nbins;
            pnegions  = fld(:,7)*nsalt*lz*area/nbins;
            pcntrions = fld(:,8)*ncntr(nvals)*lz*area/nbins;
            
            netcharge_pion = chargearr(5)*pposions;
            netcharge_nion = chargearr(6)*pnegions + cntrcharg*pcntrions;
                
            
            sumq_p = 0.5*netcharge_pion(1);
            sumq_n = 0.5*netcharge_nion(1);
           
            pos_qofr = zeros(length(netcharge_pion),1); 
            neg_qofr = zeros(length(netcharge_nion),1); 
            rqofr = zeros(length(netcharge_pion),1);
            
            pos_qofr(1,1) = sumq_p; neg_qofr(1,1) = sumq_n;
            rqofr(1,1) = 0.5*rdata(1);
            
            % Integrate charge to the edge of the brush
            for k = 2:length(netcharge_pion)
                
                sumq_p = sumq_p + 0.5*(rdata(k)-rdata(k-1))*(netcharge_pion(k)+netcharge_pion(k-1));
                sumq_n = sumq_n + 0.5*(rdata(k)-rdata(k-1))*(netcharge_nion(k)+netcharge_nion(k-1));
                pos_qofr(k,1)  = sumq_p;
                neg_qofr(k,1)  = sumq_n;
                rqofr(k,1) = 0.5*(rdata(k)+rdata(k-1));
                
            end
            
            for k = 1:length(netcharge_pion)
                
                fprintf(fnr,'%g\t%g\t%g\n',rqofr(k,1),pos_qofr(k,1),neg_qofr(k,1));
                
            end
            
            fclose(fid);
            fclose(fnr);
            
        end
        
    end
    
    % Plot Q_{z}-ion data
    
    for i = 1:nfreecnt
        
        h1 = figure;
        hold on
        box on
        set(gca,'FontSize',16)
        xlabel('$z/L_z$','FontSize',20,'Interpreter','Latex')
        ylabel('$q_{ions}(z)$','FontSize',20,'Interpreter','Latex')
        ax = zeros(4,1);
        for j = 1:4
            
            if j == 1
                dirstr = 'bl_bl';
            elseif j == 2
                dirstr = 'bl_al';
            elseif j == 3
                dirstr = 'al_bl';
            elseif j == 4
                dirstr = 'al_al';
            else
                disp('No Correct String')
                break;
            end
            
            fnr = fopen(sprintf('./All_TxtFiles_Paper/figr1_ionQz_%d_%s.txt',nfreearr(i),dirstr),'r');
            data = textscan(fnr,'%f%f%f','Headerlines',2);
            
            pl_data = zeros(length(rqofr),3);
            fld = cell2mat(data);
            pl_data(:,1) = fld(:,1);
            pl_data(:,2) = fld(:,2);
            pl_data(:,3) = fld(:,3);
            
            ax(j) = plot(pl_data(:,1)/lz,pl_data(:,2)+pl_data(:,3),'color',pclr{j},'LineWidth',2,'LineStyle',lsty{1});
%             plot(pl_data(:,1)/lz,pl_data(:,3),'color',pclr{j},'LineWidth',2,'LineStyle',lsty{2});
            fclose(fnr);
            
        end
        
        legendinfo{1} = 'Block-Block';
        legendinfo{2} = 'Block-Alter';
        legendinfo{3} = 'Alter-Block';
        legendinfo{4} = 'Alter-Alter';
        
        legend(ax,legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
        legend boxoff
        saveas(h1,sprintf('FigS_Qzion_%g',nfreearr(i)),'png');
        clear legendinfo
        
    end
    
end