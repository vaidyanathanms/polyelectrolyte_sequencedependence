%% Final Figure Plots

clear;
clc;
close all;
format long;


%% Flags

fig2a = 0 ; %fads-Npa/Npc
fig2b = 0; %delF-Npa/Npc
fig3a = 0; %delU/delN-Npa/Npc
fig3b = 0; %delS-Npa/Npc
fig4  = 0; %Density profiles @n=72
fig5  = 0; %Density profiles @n=32
old_fig4a = 0; %rhopa,rhopc-z @n=32
old_fig4b = 0; %rhopa,rhopc-z @n=100
fig6a = 1; %q(z) @n=100
fig6b = 0; %Qb vs Npa/Npc

%% Inputs

nmonfree = 30; nmongraft = 30; ngraft = 64;nsalt=510;nbackmons=10;
ncharg_chain = 15;
nfreearr = [100];
cutoff = '1.50'; lz = 120; area=53^2;
rhofree = nfreearr*nmonfree/(lz*area);


green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17]; brown=[0.6 0.2 0];
pclr = {'m',brown,green,'k','b', gold};
lsty = {'-','--',':'};
msty = {'d','s','o','x'};

nadschain = zeros(length(nfreearr),4);

%% Fig 2a. Plot adsorbed fraction as a function of number of ADSORBED CHAINS
%  Ref:adsfrac.m
if fig2a == 1
    
    fprintf('%s\n','Preparing plots for adsorbed chain fraction');
    fout = fopen(sprintf('./../../all_txtfiles/adsorbed_chain_ave_rcut_%s.dat',cutoff),'w');
    fprintf(fout,'%s\t%s\t%s\n','N_f','Arch','fraction');
    errvals = importdata('./../../err_data/errvals_fads.dat');
    
    for ncnt = 1:length(nfreearr)
        nval = nfreearr(ncnt);
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
            
            filename = sprintf('./results_adsfrac/results_%d_%s/adsfracchain_rcut_%s.lammpstrj',...
                nval,dirstr,cutoff);
            data = importdata(filename);
            nadschain(ncnt,i) = mean(data(:,3));
            fprintf(fout,'%d\t%s\t%g\n',nval,dirstr,nadschain(ncnt,i));
        end
    end
    fclose(fout);
    
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$N_{pa}/N_{pc}$','FontSize',20,'Interpreter','Latex')
    ylabel('$f_{ads}$','FontSize',20,'Interpreter','Latex')
    
    errorbar(nfreearr/ngraft,nadschain(:,1),errvals.data(:,2),'color',pclr{1},'LineWidth',2,'LineStyle',lsty{3},'Marker',msty{1},'MarkerSize',8,'MarkerFaceColor',pclr{1})
    errorbar(nfreearr/ngraft,nadschain(:,3),errvals.data(:,3),'color',pclr{3},'LineWidth',2,'LineStyle',lsty{3},'Marker',msty{3},'MarkerSize',8,'MarkerFaceColor',pclr{3})
    errorbar(nfreearr/ngraft,nadschain(:,2),errvals.data(:,4),'color',pclr{2},'LineWidth',2,'LineStyle',lsty{3},'Marker',msty{2},'MarkerSize',8,'MarkerFaceColor',pclr{2})
    errorbar(nfreearr/ngraft,nadschain(:,4),errvals.data(:,5),'color',pclr{4},'LineWidth',2,'LineStyle',lsty{3},'Marker',msty{4},'MarkerSize',8,'MarkerFaceColor',pclr{4})
    
    legendinfo{1} = 'Block-Block';
    legendinfo{2} = 'Alter-Block';
    legendinfo{3} = 'Block-Alter';
    legendinfo{4} = 'Alter-Alter';
    
    
    %overlay y = x line
    
    x = 0:0.1:1.1; y = x;
    plot(x,y,'LineWidth',2,'Color',orange,'LineStyle','--')
    
    
    legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
    legend boxoff
    saveas(h1,sprintf('./../../all_figures/adsorbchain_rcut_%s.png',cutoff));
    
end

%% Fig2b. Plot Free Energy
% Ref: whamplots.m
if fig2b == 1
    nvalsarr = [32;48;64;72;80;100];
    diff_energy = zeros(length(nvalsarr),4);
    err_energy = zeros(length(nvalsarr),4);
    fout = fopen('./../../all_txtfiles/expdeltaF_all.dat','w');
    fprintf(fout,'%s\t%s\t%s\t%s\n','Nfree','Arch','deltaF','Error');
    for nvals = 1:length(nvalsarr)
        
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
            
            fprintf('WHAM analysis for %d,\t %s\n',nvalsarr(nvals),dirstr);
            fylename = sprintf('./../../whamout_all/wham_%d_%s/whamout.txt',nvalsarr(nvals),dirstr);
            fid = fopen(fylename,'r');
            free_energy = zeros(10,3);
            header = fgetl(fid);
            k = 1;
            
            while ~feof(fid)
                
                tline = fgetl(fid);
                strarr = strsplit(tline);
                
                if strcmp(strarr{1},'#Window')
                    break;
                else
                    free_energy(k,1) = str2double(strarr{1});
                    free_energy(k,2) = str2double(strarr{2});
                    free_energy(k,3) = str2double(strarr{3});
                    k = k + 1;
                end
                
            end
            
            [minfree,indmin] = min(free_energy(:,2));
            [maxfree,indmax] = max(free_energy(:,2));
            lenfree_energy = length(free_energy(:,2));
            bulkfree_energy = mean(free_energy(floor(0.95*lenfree_energy):lenfree_energy,2));
            errbulk = 0.0; indfree = floor(0.95*lenfree_energy);
            while indfree <= lenfree_energy
                errbulk = free_energy(indfree,3)^2 + errbulk;
                indfree = indfree + 1;
            end
            
            diff_energy(nvals,i) =  minfree - bulkfree_energy;
            err_energy(nvals,i) = sqrt(free_energy(indmin,3)^2 + errbulk);
            fprintf(fout,'%d\t%s\t%g\t%g\n',nvalsarr(nvals),dirstr,diff_energy(nvals,i),err_energy(nvals,i));
            
        end
        fclose(fid);
    end
    
    
    
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$N_{pa}/N_{pc}$','FontSize',20,'Interpreter','Latex')
    ylabel('$\Delta F$ ($k_B T$)','FontSize',20,'Interpreter','Latex')
    
    for i = 1:4
        errorbar(nvalsarr/ngraft,diff_energy(:,i),err_energy(:,i),'color',pclr{i},'LineWidth',2, ...
            'LineStyle',lsty{3},'Marker',msty{i},'MarkerSize',8,'MarkerFaceColor',pclr{i})
    end
    
    legendinfo{1} = 'Block-Block';
    legendinfo{2} = 'Block-Alter';
    legendinfo{3} = 'Alter-Block';
    legendinfo{4} = 'Alter-Alter';
    
    legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
    legend boxoff
    saveas(h1,'./../../all_figures/figs_repo/delta_F','png');
    
end

%% Fig 3a. Intensive (delU/delN) with errorbar
% Ref: analyzeenergy.m

if fig3a == 1
    
    nvalsarr = [32,48,64,72];
    
    all_Emeans = fopen('./../../all_txtfiles/AllEner.dat','r');
    errUbyN = importdata('./../../all_txtfiles/deludelN.txt');
    
    varUbyN  = zeros(length(nvalsarr),4);
    meanUbyN = zeros(length(nvalsarr),4);
    
    for i = 1:3
        tline = fgetl(all_Emeans);
    end
    
    lennvals = length(nvalsarr); arrcnt = 1; archcnt = 1;
    while ~feof(all_Emeans)
        tline = fgetl(all_Emeans);
        strarr = strsplit(tline);
        meanUbyN(arrcnt,archcnt) = str2double(strarr{6});
        if rem(arrcnt,lennvals) == 0
            arrcnt = 1; archcnt = archcnt + 1;
        else
            arrcnt = arrcnt + 1;
        end
    end
    
    %Cols 3-6 (UbyN value std)
    
    kvar = 3; %See above comment
    for i = 1:4
        varUbyN(:,i) = errUbyN.data(:,kvar);
        kvar = kvar + 1;
    end
    
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$N_{pa}/N_{pc}$','FontSize',20,'Interpreter','Latex')
    ylabel('$\langle \Delta U \rangle$ ($k_B T$)','FontSize',20,'Interpreter','Latex')
    for i = 1:4
        errorbar(nvalsarr/ngraft,meanUbyN(:,i),varUbyN(:,i),'color',pclr{i},'LineWidth',2, ...
            'LineStyle',lsty{3},'Marker',msty{i},'MarkerSize',8,'MarkerFaceColor',pclr{i})
    end
    
    legendinfo{1} = 'Block-Block';
    legendinfo{2} = 'Block-Alter';
    legendinfo{3} = 'Alter-Block';
    legendinfo{4} = 'Alter-Alter';
    
    legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
    legend boxoff
    saveas(h1,'./../../all_figures/figs_repo/deltaUbyN_witherr','png');
    fclose(all_Emeans);
end

%% Fig 3b. Entropy with errorbar
% Ref: finplots_energy.m

if fig3b == 1
    nvalsarr = [32,48,64,72]; ncutoff = '1.50';
    
    fylename = './../../all_txtfiles/expdeltaF_all.dat';
    if exist(fylename, 'file') == 2
        
        all_Emeans = fopen('./../../all_txtfiles/AllEner.dat','r');
        errUbyN = importdata('./../../all_txtfiles/deludelN.txt');
        
        varUbyN  = zeros(length(nvalsarr),4);
        meanUbyN = zeros(length(nvalsarr),4);
        
        for i = 1:3
            tline = fgetl(all_Emeans);
        end
        
        lennvals = length(nvalsarr); arrcnt = 1; archcnt = 1;
        while ~feof(all_Emeans)
            tline = fgetl(all_Emeans);
            strarr = strsplit(tline);
            meanUbyN(arrcnt,archcnt) = str2double(strarr{6});
            if rem(arrcnt,lennvals) == 0
                arrcnt = 1; archcnt = archcnt + 1;
            else
                arrcnt = arrcnt + 1;
            end
        end
        
        %Cols 3-6 (N value std)
        
        kvar = 3; %See above comment
        for i = 1:4
            varUbyN(:,i) = errUbyN.data(:,kvar);
            kvar = kvar + 1;
        end
        
        
        dataN = importdata('./../../all_txtfiles/delN.txt');
        %Cols 3-6 (N value std), Cols 10-13 (N value mean)
        
        varN  = zeros(length(nvalsarr),4);
        meanN = zeros(length(nvalsarr),4);
        
        kvar = 3; kmean = 10; %See above comment
        for i = 1:4
            varN(:,i) = dataN.data(:,kvar);
            meanN(:,i) = dataN.data(:,kmean);
            kvar = kvar + 1; kmean = kmean+1;
        end
        
        
        delF = zeros(length(nvalsarr),4);
        errF = zeros(length(nvalsarr),4);
        k = 1;
        ffree = fopen(fylename,'r');
        header = fgetl(ffree);
        free_energy = zeros(length(nvalsarr)*4,2);
        
        while ~feof(ffree) && k <= length(nvalsarr)*4
            tline = fgetl(ffree);
            strarr = strsplit(tline);
            free_energy(k,1) = str2double(strarr{3});
            free_energy(k,2) = str2double(strarr{4});
            k = k + 1;
        end
        fclose(ffree);
        
        if k ~= length(nvalsarr)*4+1
            error('Unequal number of free energy columns%d\t%d\n',k,length(nvalsarr)*4+1);
        end
        
        k = 1;
        for i = 1:length(nvalsarr)
            for j = 1:4
                delF(i,j) = free_energy(k,1);
                errF(i,j) = free_energy(k,2);
                k = k + 1;
            end
        end
        
        h1 = figure;
        hold on
        box on
        set(gca,'FontSize',16)
        xlabel('$N_{pa}/N_{pc}$','FontSize',20,'Interpreter','Latex')
        ylabel('$\langle T \Delta S \rangle$ ($k_B T$)','FontSize',20,'Interpreter','Latex')
        
        for i = 1:4
            
            delS = meanUbyN(:,i) - delF(:,i);
            errfromU = varUbyN(:,i);
            errfromF = errF(:,i);
            errforS  = sqrt(errfromU.^2 + errfromF.^2);
            errorbar(nvalsarr/ngraft,delS,varUbyN(:,i),'color',pclr{i},'LineWidth',2, ...
                'LineStyle',lsty{3},'Marker',msty{i},'MarkerSize',8,'MarkerFaceColor',pclr{i})
            
            
        end
        legendinfo{1} = 'Block-Block';
        legendinfo{2} = 'Block-Alter';
        legendinfo{3} = 'Alter-Block';
        legendinfo{4} = 'Alter-Alter';
        
        legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
        legend boxoff
        saveas(h1,'deltaS','png');
        
        
    else
        fprintf('%s not found',fylename);
        error('Cannot compute Delta S');
    end
    
end

%% Fig 4. Density Profiles @n=72
% Ref: plotdens.m

if fig4 == 1
    
    nval = 72;
    fprintf('Preparing density plots for %g\n',nval);
    pclr = {'m',brown,green,'k','b', gold};
    for i = 1:4
        
        h1 = figure;
        hold on
        box on
        set(gca,'FontSize',16)
        xlabel('$r/L_z$','FontSize',20,'Interpreter','Latex')
        ylabel('$\rho(r)$','FontSize',20,'Interpreter','Latex')
        
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
        
        fprintf('Analyzing ./results_dens/results_%d_%s/dens.lammpstrj \n',nval,dirstr);
        fid = fopen(sprintf('./../../results_dens/results_%d_%s/dens.lammpstrj',nval,dirstr));
        data = textscan(fid,'%f%f%f%f%f%f%f%f','Headerlines',1);
        fld = cell2mat(data);
        rdata     = fld(:,1); %r-data
        pconnect  = fld(:,2); %neutral connector
        pneutcat  = fld(:,3); %neutral cation
        pcharcat  = fld(:,4); %charged cation
        pneutani  = fld(:,5); %neutral anion
        pcharani  = fld(:,6); %charged anion
        pposions  = fld(:,7); %positive salt/ions
        pnegions  = fld(:,8); %negative salt/ions
        
        if i == 1 || i == 3
            
            plot(rdata/lz,pconnect, 'Color', pclr{1}, 'LineWidth', 2, 'MarkerSize', 8, 'LineStyle',lsty{1});
            plot(rdata/lz,pneutcat, 'Color', pclr{2}, 'LineWidth', 2, 'MarkerSize', 8, 'LineStyle',lsty{1});
            plot(rdata/lz,pcharcat, 'Color', pclr{3}, 'LineWidth', 2, 'MarkerSize', 8, 'LineStyle',lsty{2});
            plot(rdata/lz,pneutani, 'Color', pclr{4}, 'LineWidth', 2, 'MarkerSize', 8, 'LineStyle',lsty{3});
            plot(rdata/lz,pcharani, 'Color', pclr{5}, 'LineWidth', 2, 'MarkerSize', 8, 'LineStyle',lsty{1});
            legendinfo{1} = 'Spacer';
            legendinfo{1} = 'Neutral Block PC';
            legendinfo{2} = 'Charged Block PC ';
            legendinfo{3} = 'Neutral Block PA';
            legendinfo{4} = 'Charged Block PA';
            
        else
            
            plot(rdata/lz,pconnect, 'Color', pclr{1}, 'LineWidth', 2, 'MarkerSize', 8, 'LineStyle',lsty{1});
            plot(rdata/lz,pneutcat, 'Color', pclr{2}, 'LineWidth', 2, 'MarkerSize', 8, 'LineStyle',lsty{1});
            plot(rdata/lz,pcharcat, 'Color', pclr{3}, 'LineWidth', 2, 'MarkerSize', 8, 'LineStyle',lsty{2});
            %plot(rdata/lz,pneutani, 'Color', pclr{4}, 'LineWidth', 2, 'MarkerSize', 8, 'LineStyle',lsty{1});
            plot(rdata/lz,pcharani, 'Color', pclr{5}, 'LineWidth', 2, 'MarkerSize', 8, 'LineStyle',lsty{3});
            legendinfo{1} = 'Spacer';
            legendinfo{2} = 'Neutral Block PC';
            legendinfo{3} = 'Charged Block PC ';
            %legendinfo{4} = 'Neutral Block PA';
            legendinfo{4} = 'Charged Block PA';
            
        end
        
        legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
        legend boxoff
        xlim([0 0.30])
        saveas(h1,sprintf('./../../all_figures/figs_repo/n_%d_%s',nval,dirstr),'png');
        saveas(h1,sprintf('./../../all_figures/figs_repo/n_%d_%s',nval,dirstr),'eps');
        saveas(h1,sprintf('./../../all_figures/figs_repo/n_%d_%s',nval,dirstr),'fig');
        clear legendinfo
        
    end
    
end

%% Fig 4. Density Profiles @n=32
% Ref: plotdens.m

if fig5 == 1
    
    nval = 32;
    fprintf('Preparing density plots for %g\n',nval);
    pclr = {'m',brown,green,'k','b', gold};
    for i = 1:4
        
        h1 = figure;
        hold on
        box on
        set(gca,'FontSize',16)
        xlabel('$r/L_z$','FontSize',20,'Interpreter','Latex')
        ylabel('$\rho(r)$','FontSize',20,'Interpreter','Latex')
        
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
        
        fprintf('Analyzing ./results_dens/results_%d_%s/dens.lammpstrj \n',nval,dirstr);
        fid = fopen(sprintf('./../../results_dens/results_%d_%s/dens.lammpstrj',nval,dirstr));
        data = textscan(fid,'%f%f%f%f%f%f%f%f','Headerlines',1);
        fld = cell2mat(data);

        rdata     = fld(:,1); %r-data
        pconnect  = fld(:,2); %neutral connector
        pneutcat  = fld(:,3); %neutral cation
        pcharcat  = fld(:,4); %charged cation
        pneutani  = fld(:,5); %neutral anion
        pcharani  = fld(:,6); %charged anion
        pposions  = fld(:,7); %positive salt/ions
        pnegions  = fld(:,8); %negative salt/ions

        
        if i == 1 || i == 3
            
            %plot(rdata/lz,pconnect, 'Color', pclr{1}, 'LineWidth', 2, 'MarkerSize', 8, 'LineStyle',lsty{1});
            plot(rdata/lz,pneutcat, 'Color', pclr{2}, 'LineWidth', 2, 'MarkerSize', 8, 'LineStyle',lsty{1});
            plot(rdata/lz,pcharcat, 'Color', pclr{3}, 'LineWidth', 2, 'MarkerSize', 8, 'LineStyle',lsty{2});
            plot(rdata/lz,pneutani, 'Color', pclr{4}, 'LineWidth', 2, 'MarkerSize', 8, 'LineStyle',lsty{3});
            plot(rdata/lz,pcharani, 'Color', pclr{5}, 'LineWidth', 2, 'MarkerSize', 8, 'LineStyle',lsty{1});
            %legendinfo{1} = 'Spacer';
            legendinfo{1} = 'Neutral Block PC';
            legendinfo{2} = 'Charged Block PC ';
            legendinfo{3} = 'Neutral Block PA';
            legendinfo{4} = 'Charged Block PA';
            
        else
            
            %plot(rdata/lz,pconnect, 'Color', pclr{1}, 'LineWidth', 2, 'MarkerSize', 8, 'LineStyle',lsty{1});
            plot(rdata/lz,pneutcat, 'Color', pclr{2}, 'LineWidth', 2, 'MarkerSize', 8, 'LineStyle',lsty{1});
            plot(rdata/lz,pcharcat, 'Color', pclr{3}, 'LineWidth', 2, 'MarkerSize', 8, 'LineStyle',lsty{2});
            %plot(rdata/lz,pneutani, 'Color', pclr{4}, 'LineWidth', 2, 'MarkerSize', 8, 'LineStyle',lsty{1});
            plot(rdata/lz,pcharani, 'Color', pclr{5}, 'LineWidth', 2, 'MarkerSize', 8, 'LineStyle',lsty{3});
            %legendinfo{1} = 'Spacer';
            legendinfo{1} = 'Neutral Block PC';
            legendinfo{2} = 'Charged Block PC ';
            %legendinfo{3} = 'Neutral Block PA';
            legendinfo{3} = 'Charged Block PA';
            
        end
        
        legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
        legend boxoff
        xlim([0 0.30])
        saveas(h1,sprintf('./../../all_figures/figs_repo/n_%d_%s',nval,dirstr),'png');
        saveas(h1,sprintf('./../../all_figures/figs_repo/n_%d_%s',nval,dirstr),'eps');
        saveas(h1,sprintf('./../../all_figures/figs_repo/n_%d_%s',nval,dirstr),'fig');
        clear legendinfo
        
    end
    
end

%% Fig old_4a. Density Profiles @n-32
% Ref:paper_methods.m

if old_fig4a == 1
    
    nfree = 32;
    ncntr = abs(ngraft*nmongraft/2-nmonfree*nfree/2);
    ntotmons = nfree*nmonfree+ngraft*(nmongraft+nbackmons) + ncntr+nsalt*2;
    
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
        
        fid = fopen(sprintf('./../../results_dens/results_%d_%s/grpdens.lammpstrj',nfree,dirstr));
        data = textscan(fid,'%f%f%f','Headerlines',1);
        fld = cell2mat(data);
        
        rdata     = fld(:,1);
        pl_data(:,1) = fld(:,2)*ngraft*nmongraft/ntotmons;
        pl_data(:,2) = fld(:,3)*nfree*nmonfree/ntotmons;
        
        patch(rdata/lz,pl_data(:,1),gold);
        alpha(0.3);
        patch(rdata/lz,pl_data(:,2),'b')
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
    legendinfo{1} = '$\rho_{pc}(z)$';
    legendinfo{2} = '$\rho_{pa}(z)$';
    legend(h(1), legendinfo, 'Interpreter','Latex','FontSize',16)
    box(legend,'off')
end

%% Fig old_4b. Density Profiles
% Ref:paper_methods.m

if old_fig4b == 1
    
    nfree = 100;
    ncntr = abs(ngraft*nmongraft/2-nmonfree*nfree/2);
    ntotmons = nfree*nmonfree+ngraft*(nmongraft+nbackmons) + ncntr+nsalt*2;
    
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
        
        fid = fopen(sprintf('./../../results_dens/results_%d_%s/grpdens.lammpstrj',nfree,dirstr));
        data = textscan(fid,'%f%f%f','Headerlines',1);
        fld = cell2mat(data);
        
        rdata     = fld(:,1);
        pl_data(:,1) = fld(:,2)*ngraft*nmongraft/ntotmons;
        pl_data(:,2) = fld(:,3)*nfree*nmonfree/ntotmons;
        
        patch(rdata/lz,pl_data(:,1),gold);
        alpha(0.3);
        patch(rdata/lz,pl_data(:,2),'b')
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
    legendinfo{1} = '$\rho_{pc}(z)$';
    legendinfo{2} = '$\rho_{pa}(z)$';
    legend(h(1), legendinfo, 'Interpreter','Latex','FontSize',16)
    box(legend,'off')
    
end
%% Fig 6a. Charge density, q(z)
% Ref:densananew.m

if fig6a == 1
    
    fprintf('%s\n','Preparing plots for charge density');
    nspacer = 10;
    nfree = 80;rbin = 1;nbins=lz/rbin;
    ntotal = (nfreearr + ngraft)*nmonfree + 2*nsalt + nbackmons*ngraft;
    chargearr = [0;0;1;0;-1;1;-1];
    
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$z/L_z$','FontSize',20,'Interpreter','Latex')
    ylabel('$q(z)$','FontSize',20,'Interpreter','Latex')
    
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

        nposions = nsalt + nfree*ncharg_chain;
        nnegions = nsalt + ngraft*ncharg_chain;
        
        fprintf('Preparing plots for %s\n',dirstr);
        fid = fopen(sprintf('./../../results_dens/results_%d_%s/dens.lammpstrj',nfree,dirstr));
        data = textscan(fid,'%f%f%f%f%f%f%f%f','Headerlines',1);
        fld = cell2mat(data);
        
        rdata     = fld(:,1);
        pbackbon_brush = fld(:,2)*nbackmons*ngraft*lz*area/nbins;
        pneutral_brush = fld(:,3)*ngraft*nmongraft/2*lz*area/nbins;
        pnegativ_brush = fld(:,4)*ngraft*nmongraft/2*lz*area/nbins;
        pneutral_free  = fld(:,5)*nfree*nmonfree/2*lz*area/nbins;
        ppositive_free = fld(:,6)*nfree*nmonfree/2*lz*area/nbins;
        pposions  = fld(:,7)*nposions*lz*area/nbins;
        pnegions  = fld(:,8)*nnegions*lz*area/nbins;
        
        
        
        
        netcharge = chargearr(1)*pbackbon_brush + chargearr(2)*pneutral_brush ...
            + chargearr(3)*pnegativ_brush + chargearr(4)*pneutral_free + ...
            chargearr(5)*ppositive_free + chargearr(6)*pposions + chargearr(7)*pnegions;
        
        %trapz(rdata,pposions)
        
        
        
        intdata = zeros(length(rdata),1);
        sumval  = 0.5*netcharge(1); 
        intdata(1) = sumval;
        for icnt = 2:length(rdata)
            intdata(icnt) = intdata(icnt-1) + 0.5*(rdata(icnt)-rdata(icnt-1))*(netcharge(icnt) + netcharge(icnt-1));
        end
          
        
        plot(rdata/lz,netcharge, 'color', pclr{i}, 'LineWidth', 2)
        
    end
    legendinfo{1} = 'Block-Block';
    legendinfo{2} = 'Block-Alter';
    legendinfo{3} = 'Alter-Block';
    legendinfo{4} = 'Alter-Alter';
    
    legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
    legend boxoff
    saveas(h1,sprintf('./../../all_figures/n%d_integral_netcharge',nfree),'png');
    
end

%% Fig 6b: Qb - Npc/Npa

if fig6b == 1
    
    fprintf('%s\n','Preparing net bound charge plot');
    chargearr = [0;0;1;0;-1;1;-1];
    rbin = 1;nbins=lz/rbin;
    ncntr = abs((nfreearr*nmonfree-ngraft*nmongraft)/2);
    ntotarr = (nfreearr + ngraft)*nmonfree + 2*nsalt + nbackmons*ngraft;

    
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$r$','FontSize',20,'Interpreter','Latex')
    ylabel('$f$','FontSize',20,'Interpreter','Latex')
    for nvals = 1:length(nfreearr)  
        
        nfree = nfreearr(nvals);
        
        fnr = fopen(sprintf('./../../all_txtfiles/fig6b_QnetBound_%d.txt',nfreearr(nvals)),'w');
        fprintf(fnr,'%s \n','NetCharge: Q_{b}=\Delta(n_g)*\int(\sum(q_j n_j(z)dz, j=all entities)z=0,Lz))');
        fprintf(fnr,'%s\t %s\n','Arch','Q_{b}');
        
        
       
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
            
            nposions = nsalt + nfree*ncharg_chain;
            nnegions = nsalt + ngraft*ncharg_chain;
            
            fprintf('%s\t %d\t %d\t %d\n',dirstr,nfree,nposions,nnegions);

            
            fid = fopen(sprintf('./../../results_dens/results_%d_%s/dens.lammpstrj',nfree,dirstr));
            if fid <= 0
                fprintf('Did not find dens.lammpstrj for %d\t%s\n', nfree,dirstr);
            end
            data = textscan(fid,'%f%f%f%f%f%f%f%f','Headerlines',1);
            fld = cell2mat(data);
       
            rdata = fld(:,1);
            pbackbon_brush = fld(:,2)*ngraft*nbackmons*lz*area/nbins;
            pneutral_brush = fld(:,3)*ngraft*nmongraft/2*lz*area/nbins;
            pnegativ_brush = fld(:,4)*ngraft*nmongraft/2*lz*area/nbins;
            pneutral_free  = fld(:,5)*nfree*nmonfree/2*lz*area/nbins;
            ppositive_free = fld(:,6)*nfree*nmonfree/2*lz*area/nbins;
            pposions  = fld(:,7)*nposions*lz*area/nbins;
            pnegions  = fld(:,8)*nnegions*lz*area/nbins;
            
            
            netcharge = chargearr(1)*pbackbon_brush + chargearr(2)*pneutral_brush ...
                + chargearr(3)*pnegativ_brush + chargearr(4)*pneutral_free + ...
                chargearr(5)*ppositive_free + chargearr(6)*pposions + chargearr(7)*pnegions;
            
            sumq = 0.5*netcharge(1);
            
            fid_g  = fopen(sprintf('./../../results_dens/results_%d_%s/grpdens.lammpstrj',nfree,dirstr));
            data_g = textscan(fid_g,'%f%f%f','Headerlines',1);
            
            fld_g   = cell2mat(data_g);
            dens_g  = fld_g(:,2);
            [maxden, imaxden]  = max(dens_g);
            
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
            %fprintf('%s\t%g\n',dirstr,sumq);
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
        
        fnr = fopen(sprintf('./../../all_txtfiles/fig6b_QnetBound_%d.txt',nfreearr(i)),'r');
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
    saveas(h1,'./../../all_figures/Fig1b_QnetBound_%s','png');
    clear legendinfo
    
end