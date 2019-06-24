%% To compute the adsorbed values in different formats

clear all
close all
clc
format long


%% Flags

mrdf = 0;
m3 = 0;
m4 = 0;
m5 = 0;
m6 = 0;
mads = 1;

%% Inputs

nmonfree = 30; nmongraft = 30; ngraft = 64;
nfreearr = [16;32;48;64;80;100;150];
cutoff = 0.9; lz = 120; area=53^2;
rhofree = nfreearr*30/(lz*area);

green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17];
pclr = {'r','b',green,'k','m', gold};
lsty = {'-','--',':'};
msty = {'d','s','o','x'};

%% Method 1: int(4*pi*r^2*g(r)dr,0,rcut=1.5sigma)/ng

rcut = 1.5;

if mrdf == 1
    for nvals = 1:length(nfreearr)
        
        nfree = nfreearr(nvals);
        fnr = fopen(sprintf('nofr_%d_%g.txt',nfreearr(nvals),rcut),'w');
        fprintf(fnr,'%s \n','Method 1: int(4*pi*r^2*g(r)dr,0,rcut=1.5sigma)');
        fprintf(fnr,'%s\t %s\n','Arch','n(r_cut)');
        
        h1 = figure;
        hold on
        box on
        set(gca,'FontSize',16)
        xlabel('$r$','FontSize',20,'Interpreter','Latex')
        ylabel('$n_f(r)$','FontSize',20,'Interpreter','Latex')
        
        for i = 1:4
            
            if i == 1
                dirstr = 'blbl';
            elseif i == 2
                dirstr = 'blal';
            elseif i == 3
                dirstr = 'albl';
            else
                dirstr = 'alal';
            end
            
            fid = fopen(sprintf('./n_%d/rdf_%s.txt',nfree,dirstr));
            data = textscan(fid,'%f%f','Headerlines',1);
            fld = cell2mat(data);
            rdata    = fld(:,1);
            rdfdata  = fld(:,2);
            nofrdata = zeros(length(rdata)-1,1);
            rnrdata = zeros(length(rdata)-1,1); 
            rsqgrvals = (rdata.^2).*(rdfdata);
            data1 = rdata(1:32); data2=rsqgrvals(1:32);
%             trapz(data1,data2)*4*pi*nmonfree*rhofree(nvals)
            nofrdata(1,1) = rsqgrvals(1,1)*0.5;
            rnrdata(1,1) = rdata(1,1)*0.5;
            for j = 2:length(rdata)-1
                nofrdata(j,1) = nofrdata(j-1,1) + 0.5*(rdata(j)-rdata(j-1))*(rsqgrvals(j-1)+rsqgrvals(j));
                rnrdata(j,1)  = 0.5*(rdata(j)+rdata(j-1));
            end
            
            plot(rnrdata,4*pi*rhofree(nvals)*nofrdata,'color',pclr{i},'LineStyle',lsty{nvals},'LineWidth',2)
            
            rspline = 0:0.01:max(rnrdata);
            nofrspline = spline(rnrdata,nofrdata,rspline);
            nofrval = 4*pi*rhofree(nvals)*spline(rspline,nofrspline,rcut);
            fprintf(fnr,'%s\t%g\n',dirstr,nofrval);
            
            fclose(fid);
            
        end
        
        legendinfo{1} = 'Block-Block';
        legendinfo{2} = 'Block-Alter';
        legendinfo{3} = 'Alter-Block';
        legendinfo{4} = 'Alter-Alter';
        
        legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
        legend boxoff
        xlim([0 2.5] )
        saveas(h1,sprintf('PEnofr_%d',nfree),'png');
        fclose(fnr);
        
    end
    
    clear legendinfo
    
    % Plot the n(r) values
    
    data_bb = zeros(length(nfreearr),1);
    data_ab = zeros(length(nfreearr),1);
    data_ba = zeros(length(nfreearr),1);
    data_aa = zeros(length(nfreearr),1);
    
    for i = 1:length(nfreearr)
        
        fnr = fopen(sprintf('nofr_%d_%g.txt',nfreearr(i),rcut));
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
    saveas(h1,'Method1','png');
end


%% Method 2: int(4*pi*r^2*g(r)dr,0,rcut=2.0sigma)/ng

rcut = 2.0;

if mrdf == 1
    for nvals = 1:length(nfreearr)
        
        nfree = nfreearr(nvals);
        fnr = fopen(sprintf('nofr_%d_%g.txt',nfreearr(nvals),rcut),'w');
        fprintf(fnr,'%s \n','Method 2: int(4*pi*r^2*g(r)dr,0,rcut=2.0sigma)');
        fprintf(fnr,'%s\t %s\n','Arch','n(r_cut)');
        
        h1 = figure;
        hold on
        box on
        set(gca,'FontSize',16)
        xlabel('$r$','FontSize',20,'Interpreter','Latex')
        ylabel('$n_f(r)$','FontSize',20,'Interpreter','Latex')
        
        for i = 1:4
            
            if i == 1
                dirstr = 'blbl';
            elseif i == 2
                dirstr = 'blal';
            elseif i == 3
                dirstr = 'albl';
            else
                dirstr = 'alal';
            end
            
            fid = fopen(sprintf('./n_%d/rdf_%s.txt',nfree,dirstr));
            data = textscan(fid,'%f%f','Headerlines',1);
            fld = cell2mat(data);
            rdata    = fld(:,1);
            rdfdata  = fld(:,2);
            nofrdata = zeros(length(rdata)-1,1);
            rnrdata = zeros(length(rdata)-1,1); 
            rsqgrvals = (rdata.^2).*(rdfdata);
            data1 = rdata(1:32); data2=rsqgrvals(1:32);
%             trapz(data1,data2)*4*pi*nmonfree*rhofree(nvals)
            nofrdata(1,1) = rsqgrvals(1,1)*0.5;
            rnrdata(1,1) = rdata(1,1)*0.5;
            for j = 2:length(rdata)-1
                nofrdata(j,1) = nofrdata(j-1,1) + 0.5*(rdata(j)-rdata(j-1))*(rsqgrvals(j-1)+rsqgrvals(j));
                rnrdata(j,1)  = 0.5*(rdata(j)+rdata(j-1));
            end
            
            plot(rnrdata,4*pi*rhofree(nvals)*nofrdata,'color',pclr{i},'LineStyle',lsty{nvals},'LineWidth',2)
            
            rspline = 0:0.01:max(rnrdata);
            nofrspline = spline(rnrdata,nofrdata,rspline);
            nofrval = 4*pi*rhofree(nvals)*spline(rspline,nofrspline,rcut);
            fprintf(fnr,'%s\t%g\n',dirstr,nofrval);
            
            fclose(fid);
            
        end
        
        legendinfo{1} = 'Block-Block';
        legendinfo{2} = 'Block-Alter';
        legendinfo{3} = 'Alter-Block';
        legendinfo{4} = 'Alter-Alter';
        
        legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
        legend boxoff
        xlim([0 2.5] )
        saveas(h1,sprintf('PEnofr_%d',nfree),'png');
        fclose(fnr);
        
    end
    
    clear legendinfo
    
    % Plot the n(r) values
    
    data_bb = zeros(length(nfreearr),1);
    data_ab = zeros(length(nfreearr),1);
    data_ba = zeros(length(nfreearr),1);
    data_aa = zeros(length(nfreearr),1);
    
    for i = 1:length(nfreearr)
        
        fnr = fopen(sprintf('nofr_%d_%g.txt',nfreearr(i),rcut));
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
    saveas(h1,'Method2','png');
end
    

%% Method 3: int(rhop(z)dz,0,0.9rhogmax)/ng

if m3 == 1
    for nvals = 1:length(nfreearr)
        
        nfree = nfreearr(nvals);
        fnr = fopen(sprintf('fMethod3_%d.txt',nfreearr(nvals)),'w');
        fprintf(fnr,'%s \n','Method 3: int(rhop(z)dz,0,0.9rhogmax)/ng');
        fprintf(fnr,'Arch\t ht (cutoff = %g)\t  FreeIntegral\t Rat_Wrt_Ng\n',cutoff);
        
        h4 = figure;
        hold on
        box on
        set(gca,'FontSize',16)
        xlabel('$r/L_z$','FontSize',20,'Interpreter','Latex')
        ylabel('$\rho^{g}(r)$','FontSize',20,'Interpreter','Latex')
        
        for i = 1:4
            
            if i == 1
                dirstr = 'blbl';
            elseif i == 2
                dirstr = 'blal';
            elseif i == 3
                dirstr = 'albl';
            else
                dirstr = 'alal';
            end
            
            fid = fopen(sprintf('./n_%d/grp_%s.txt',nfree,dirstr));
            data = textscan(fid,'%f%f%f','Headerlines',1);
            fld = cell2mat(data);
            zdata  = fld(:,1);
            pepos  = fld(:,2);
            peneg  = fld(:,3);
            plot(zdata/lz,pepos,'color',pclr{i},'LineWidth',2,'LineStyle',lsty{nvals})
            
            maxdenval = max(pepos); cutoffval = (1-cutoff)*maxdenval;
            
            %spline fit
            zspline = 0:0.01:max(zdata);
            denspline = spline(zdata,pepos,zspline);
            pval = 0;
            for j = 1:length(zspline)-1
                if(denspline(j+1) <= cutoffval && denspline(j) >= cutoffval)
                    pval = j;
                    break;
                end
            end
            if pval == 0
                disp('Could not find the right height')
            end
            
            %nchk = (53^2)*(7120/1920)*trapz(zdata,pegraft)/1920
            brushspline_pos = spline(zdata,pepos,zspline);
            freespline_neg = spline(zdata,peneg,zspline);
            integpos = 0.0;
            integneg = 0.0;
            for j = 2:pval
                integpos = integpos+0.5*(zspline(j)-zspline(j-1))*(brushspline_pos(j)+brushspline_pos(j+1));
                integneg = integneg+0.5*(zspline(j)-zspline(j-1))*(freespline_neg(j)+freespline_neg(j+1));
            end
            
            integpos = area*integpos*ngraft*15;
            integneg = area*integneg*nfree*15;
            rat = integneg/(ngraft*nmongraft);
            ht_cut = 0.5*(zspline(pval)+zspline(pval+1));
            fprintf(fnr,'%s\t%g\t%g\t%g\n',dirstr, ht_cut, integneg, rat);
            fclose(fid);
            
            
        end
        
        legendinfo{1} = 'Block-Block';
        legendinfo{2} = 'Block-Alter';
        legendinfo{3} = 'Alter-Block';
        legendinfo{4} = 'Alter-Alter';
        
        legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
        legend boxoff
        saveas(h4,sprintf('Method3_%d',nfree),'png');
        
        fclose(fnr);
        
    end
    
    
    % Plot the f values
    
    data_bb = zeros(length(nfreearr),1);
    data_ab = zeros(length(nfreearr),1);
    data_ba = zeros(length(nfreearr),1);
    data_aa = zeros(length(nfreearr),1);
    
    for i = 1:length(nfreearr)
        
        fnr = fopen(sprintf('fMethod3_%d.txt',nfreearr(i)));
        data = textscan(fnr,'%s%f%f%f','Headerlines',2);
        
        fld = cell2mat(data(4));
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
    saveas(h1,'Method3_Integral_Free_Over_Ng','png');
end

%% Method 4: int(rhop-(z)dz,0,0.9rhogmax)/ng

if m4 == 1
    for nvals = 1:length(nfreearr)
        
        nfree = nfreearr(nvals);
        fnr = fopen(sprintf('fMethod4_%d.txt',nfreearr(nvals)),'w');
        fprintf(fnr,'%s \n','Method4: int(rhop-(z)dz,0,0.9rhogmax)/ng');
        fprintf(fnr,'Arch\t ht (cutoff = %g)\t Brush(PosCharge)Integral\t RatWrtNg\n',cutoff);
        
        h4 = figure;
        hold on
        box on
        set(gca,'FontSize',16)
        xlabel('$r/L_z$','FontSize',20,'Interpreter','Latex')
        ylabel('$\rho^{g}(r)$','FontSize',20,'Interpreter','Latex')
        
        for i = 1:4
            
            if i == 1
                dirstr = 'blbl';
            elseif i == 2
                dirstr = 'blal';
            elseif i == 3
                dirstr = 'albl';
            else
                dirstr = 'alal';
            end
            
            fid = fopen(sprintf('./n_%d/grp_%s.txt',nfree,dirstr));
            data = textscan(fid,'%f%f%f','Headerlines',1);
            fld = cell2mat(data);
            zdata  = fld(:,1);
            pepos  = fld(:,2);
            fid2 = fopen(sprintf('./n_%d/dens_%s.txt',nfree,dirstr));
            data2 = textscan(fid2,'%f%f%f','Headerlines',1);
            peneg  = fld(:,3);
            plot(zdata/lz,pepos,'color',pclr{i},'LineWidth',2,'LineStyle',lsty{nvals})
            
            maxdenval = max(pepos); cutoffval = (1-cutoff)*maxdenval;
            
            %spline fit
            zspline = 0:0.01:max(zdata);
            denspline = spline(zdata,pepos,zspline);
            pval = 0;
            for j = 1:length(zspline)-1
                if(denspline(j+1) <= cutoffval && denspline(j) >= cutoffval)
                    pval = j;
                    break;
                end
            end
            if pval == 0
                disp('Could not find the right height')
            end
            
            %nchk = (53^2)*(7120/1920)*trapz(zdata,pegraft)/1920
            brushspline_pos = spline(zdata,pepos,zspline);
            freespline_neg = spline(zdata,peneg,zspline);
            integpos = 0.0;
            integneg = 0.0;
            for j = 2:pval
                integpos = integpos+0.5*(zspline(j)-zspline(j-1))*(brushspline_pos(j)+brushspline_pos(j+1));
                integneg = integneg+0.5*(zspline(j)-zspline(j-1))*(freespline_neg(j)+freespline_neg(j+1));
            end
            
            integpos = area*integpos*ngraft*15;
            integneg = area*integneg*nfree*15;
            rat = integneg/(nmongraft*ngraft);
            ht_cut = 0.5*(zspline(pval)+zspline(pval+1));
            fprintf(fnr,'%s\t%g\t%g\t%g\n',dirstr, ht_cut, integpos, rat);
            fclose(fid);
            
            
        end
        
        legendinfo{1} = 'Block-Block';
        legendinfo{2} = 'Block-Alter';
        legendinfo{3} = 'Alter-Block';
        legendinfo{4} = 'Alter-Alter';
        
        legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
        legend boxoff
        saveas(h4,sprintf('Method4_%d',nfree),'png');
        
        fclose(fnr);
        
    end
    
    
    % Plot the f values
    
    data_bb = zeros(length(nfreearr),1);
    data_ab = zeros(length(nfreearr),1);
    data_ba = zeros(length(nfreearr),1);
    data_aa = zeros(length(nfreearr),1);
    
    for i = 1:length(nfreearr)
        
        fnr = fopen(sprintf('fMethod4_%d.txt',nfreearr(i)));
        data = textscan(fnr,'%s%f%f%f','Headerlines',2);
        
        fld = cell2mat(data(4));
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
    saveas(h1,'Method4_Integral_rhopminus','png');
end

%% Method 5: int(rhop(z)dz,0,0.9rhogmax)/int(rhog(z)dz,0,0.9rhogmax)

if m5 == 1
    for nvals = 1:length(nfreearr)
        
        nfree = nfreearr(nvals);
        fnr = fopen(sprintf('fMethod5_%d.txt',nfreearr(nvals)),'w');
        fprintf(fnr,'%s \n','Method5: int(rhop(z)dz,0,0.9rhogmax)/int(rhog(z)dz,0,0.9rhogmax)');
        fprintf(fnr,'Arch\t ht (cutoff = %g)\t Brush(All)Integral\t Free(All)Integral\t Rat\n',cutoff);
        
        h4 = figure;
        hold on
        box on
        set(gca,'FontSize',16)
        xlabel('$r/L_z$','FontSize',20,'Interpreter','Latex')
        ylabel('$\rho^{g}(r)$','FontSize',20,'Interpreter','Latex')
        
        for i = 1:4
            
            if i == 1
                dirstr = 'blbl';
            elseif i == 2
                dirstr = 'blal';
            elseif i == 3
                dirstr = 'albl';
            else
                dirstr = 'alal';
            end
            
            fid = fopen(sprintf('./n_%d/grp_%s.txt',nfree,dirstr));
            data = textscan(fid,'%f%f%f','Headerlines',1);
            fld = cell2mat(data);
            zdata  = fld(:,1);
            pepos  = fld(:,2);
            peneg  = fld(:,3);
            plot(zdata/lz,pepos,'color',pclr{i},'LineWidth',2,'LineStyle',lsty{nvals})
            
            maxdenval = max(pepos); cutoffval = (1-cutoff)*maxdenval;
            
            %spline fit
            zspline = 0:0.01:max(zdata);
            denspline = spline(zdata,pepos,zspline);
            pval = 0;
            for j = 1:length(zspline)-1
                if(denspline(j+1) <= cutoffval && denspline(j) >= cutoffval)
                    pval = j;
                    break;
                end
            end
            if pval == 0
                disp('Could not find the right height')
            end
            
            %nchk = (53^2)*(7120/1920)*trapz(zdata,pegraft)/1920
            brushspline_pos = spline(zdata,pepos,zspline);
            freespline_neg = spline(zdata,peneg,zspline);
            integpos = 0.0;
            integneg = 0.0;
            for j = 2:pval
                integpos = integpos+0.5*(zspline(j)-zspline(j-1))*(brushspline_pos(j)+brushspline_pos(j+1));
                integneg = integneg+0.5*(zspline(j)-zspline(j-1))*(freespline_neg(j)+freespline_neg(j+1));
            end
            
            integpos = area*integpos*ngraft*15;
            integneg = area*integneg*nfree*15;
            rat = integneg/integpos;
            ht_cut = 0.5*(zspline(pval)+zspline(pval+1));
            fprintf(fnr,'%s\t%g\t%g\t%g\t%g\n',dirstr, ht_cut, integpos, integneg, rat);
            fclose(fid);
            
            
        end
        
        legendinfo{1} = 'Block-Block';
        legendinfo{2} = 'Block-Alter';
        legendinfo{3} = 'Alter-Block';
        legendinfo{4} = 'Alter-Alter';
        
        legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
        legend boxoff
        saveas(h4,sprintf('Method5_%d',nfree),'png');
        
        fclose(fnr);
        
    end
    
    
    % Plot the f values
    
    data_bb = zeros(length(nfreearr),1);
    data_ab = zeros(length(nfreearr),1);
    data_ba = zeros(length(nfreearr),1);
    data_aa = zeros(length(nfreearr),1);
    
    for i = 1:length(nfreearr)
        
        fnr = fopen(sprintf('fMethod5_%d.txt',nfreearr(i)));
        data = textscan(fnr,'%s%f%f%f%f','Headerlines',2);
        
        fld = cell2mat(data(5));
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
    saveas(h1,'Method5_Rat_Freeall_Brushall_WithCutoff','png');
end

%% Method 6: int(rhop-(z)dz,0,0.9rhogmax)/int(rhog+(z)dz,0,0.9rhogmax)

if m6 == 1
    for nvals = 1:length(nfreearr)
        
        nfree = nfreearr(nvals);
        fnr = fopen(sprintf('fMikeMethod_%d.txt',nfreearr(nvals)),'w');
        fprintf(fnr,'%s \n','Method6: int(rhop-(z)dz,0,0.9rhogmax)/int(rhog+(z)dz,0,0.9rhogmax)');
        fprintf(fnr,'Arch\t ht (cutoff = %g)\t Brush(PosCharge)Integral\t Free(NegCharge)Integral\t Rat\n',cutoff);
        
        h4 = figure;
        hold on
        box on
        set(gca,'FontSize',16)
        xlabel('$r/L_z$','FontSize',20,'Interpreter','Latex')
        ylabel('$\rho^{g}(r)$','FontSize',20,'Interpreter','Latex')
        
        for i = 1:4
            
            if i == 1
                dirstr = 'blbl';
            elseif i == 2
                dirstr = 'blal';
            elseif i == 3
                dirstr = 'albl';
            else
                dirstr = 'alal';
            end
            
            fid = fopen(sprintf('./n_%d/dens_%s.txt',nfree,dirstr));
            data = textscan(fid,'%f%f%f','Headerlines',1);
            fld = cell2mat(data);
            zdata  = fld(:,1);
            pepos  = fld(:,2);
            peneg  = fld(:,3);
            plot(zdata/lz,pepos,'color',pclr{i},'LineWidth',2,'LineStyle',lsty{nvals})
            
            maxdenval = max(pepos); cutoffval = (1-cutoff)*maxdenval;
            
            %spline fit
            zspline = 0:0.01:max(zdata);
            denspline = spline(zdata,pepos,zspline);
            pval = 0;
            for j = 1:length(zspline)-1
                if(denspline(j+1) <= cutoffval && denspline(j) >= cutoffval)
                    pval = j;
                    break;
                end
            end
            if pval == 0
                disp('Could not find the right height')
            end
            
            %nchk = (53^2)*(7120/1920)*trapz(zdata,pegraft)/1920
            brushspline_pos = spline(zdata,pepos,zspline);
            freespline_neg = spline(zdata,peneg,zspline);
            integpos = 0.0;
            integneg = 0.0;
            for j = 2:pval
                integpos = integpos+0.5*(zspline(j)-zspline(j-1))*(brushspline_pos(j)+brushspline_pos(j+1));
                integneg = integneg+0.5*(zspline(j)-zspline(j-1))*(freespline_neg(j)+freespline_neg(j+1));
            end
            
            integpos = area*integpos*ngraft*15;
            integneg = area*integneg*nfree*15;
            rat = integneg/integpos;
            ht_cut = 0.5*(zspline(pval)+zspline(pval+1));
            fprintf(fnr,'%s\t%g\t%g\t%g\t%g\n',dirstr, ht_cut, integpos, integneg, rat);
            fclose(fid);
            
            
        end
        
        legendinfo{1} = 'Block-Block';
        legendinfo{2} = 'Block-Alter';
        legendinfo{3} = 'Alter-Block';
        legendinfo{4} = 'Alter-Alter';
        
        legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
        legend boxoff
        saveas(h4,sprintf('mikemethod_%d',nfree),'png');
        
        fclose(fnr);
        
    end
    
    
    % Plot the f values
    
    data_bb = zeros(length(nfreearr),1);
    data_ab = zeros(length(nfreearr),1);
    data_ba = zeros(length(nfreearr),1);
    data_aa = zeros(length(nfreearr),1);
    
    for i = 1:length(nfreearr)
        
        fnr = fopen(sprintf('fMikeMethod_%d.txt',nfreearr(i)));
        data = textscan(fnr,'%s%f%f%f%f','Headerlines',2);
        
        fld = cell2mat(data(5));
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
    saveas(h1,'Method6_Mike','png');
end


%% Method 7: 1/n_h \sum(all free such that it is within cutoff of at least 1 brush)

str_rcutarr = {'1.10';'1.50';'2.00'};
rcutarr = [1.10;1.50;2.00];

if mads == 1
    
    for rcutvals = 1:length(str_rcutarr)
        
        for nvals = 1:length(nfreearr)
        
            nfree = nfreearr(nvals);
            fnr = fopen(sprintf('n_%d/adsfracavg_%s.txt',nfreearr(nvals),str_rcutarr{rcutvals}),'w');
            fprintf(fnr,'Method 7: Sum of all free mons which is within cutoff of atleast of 1 brush, rcut= %s\n',str_rcutarr{rcutvals});
            fprintf(fnr,'%s\t %s\n','Arch','f(r_ads)');
        
            for i = 1:4
            
                if i == 1
                    dirstr = 'blbl';
                elseif i == 2
                    dirstr = 'blal';
                elseif i == 3
                    dirstr = 'albl';
                else
                    dirstr = 'alal';
                end

                
                data = importdata(sprintf('./n_%d/n%d_results/adsfrac_rcut_%s_%s.txt',nfree,nfree,str_rcutarr{rcutvals},dirstr));
                fld  = data;
                rdata  = fld(:,1);
                fdata  = fld(:,2);
                fmean  = mean(fdata);
                fprintf(fnr,'%s\t%g\n',dirstr,fmean);
        
            end
        
            fclose(fnr);
        
        end
       
    
        % Plot the f_{ads} values
    
        data_bb = zeros(length(nfreearr),1);
        data_ab = zeros(length(nfreearr),1);
        data_ba = zeros(length(nfreearr),1);
        data_aa = zeros(length(nfreearr),1);
    
        for i = 1:length(nfreearr)
        
            fnr  = fopen(sprintf('n_%d/adsfracavg_%s.txt',nfreearr(i),str_rcutarr{rcutvals}));
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
        ylabel('$f_{ads}$','FontSize',20,'Interpreter','Latex')
        
        plot(nfreearr/ngraft,data_bb/(nmongraft*ngraft),'color',pclr{1},'LineWidth',2,'LineStyle',lsty{3},'Marker',msty{1},'MarkerSize',8,'MarkerFaceColor',pclr{1})
        plot(nfreearr/ngraft,data_ba/(nmongraft*ngraft),'color',pclr{2},'LineWidth',2,'LineStyle',lsty{3},'Marker',msty{2},'MarkerSize',8,'MarkerFaceColor',pclr{2})
        plot(nfreearr/ngraft,data_ab/(nmongraft*ngraft),'color',pclr{3},'LineWidth',2,'LineStyle',lsty{3},'Marker',msty{3},'MarkerSize',8,'MarkerFaceColor',pclr{3})
        plot(nfreearr/ngraft,data_aa/(nmongraft*ngraft),'color',pclr{4},'LineWidth',2,'LineStyle',lsty{3},'Marker',msty{4},'MarkerSize',8,'MarkerFaceColor',pclr{4})
        
    
        legendinfo{1} = 'Block-Block';
        legendinfo{2} = 'Block-Alter';
        legendinfo{3} = 'Alter-Block';
        legendinfo{4} = 'Alter-Alter';
        
        legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
        legend boxoff
        saveas(h1,sprintf('Method7_rcut_%s',str_rcutarr{rcutvals}),'png');
        
    end
    
end