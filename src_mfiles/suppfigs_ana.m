%% Supplementary Information Figures plot

clear;
clc;
close all;
format long;

%% Flags

rdfplot = 0;
delUplot = 1;
delFzplot = 0;
histplot = 0;
iondens = 0;


%% Inputs

nmonfree = 30; nmongraft = 30; ngraft = 64;nsalt=510;nbackmons=10;
nbase = nbackmons*ngraft;
nfreearr = [16;32;48;64;80;100;150];
cutoff = '1.50'; lz = 120; area=53^2;
rhofree = nfreearr*nmonfree/(lz*area);


green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17]; brown=[0.6 0.2 0];
pclr = {'m',brown,green,'k','b', gold};
lsty = {'-','--',':'};
msty = {'d','s','o','x'};

nadschain = zeros(length(nfreearr),4);

%% Plot RDF for n=64

if rdfplot == 1
    nfree = 64;
    for j = 1:length(nfree)
        nval = nfree(j);
        h1 = figure;
        hold on
        box on
        set(gca,'FontSize',16)
        xlabel('$r$ ($\sigma$)','FontSize',20,'Interpreter','Latex')
        ylabel('$g_{pa-pc}$($r$)','FontSize',20,'Interpreter','Latex')
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
            fid = fopen(sprintf('./n_%d/results_%d_%s/PErdf.lammpstrj',nval,nval,dirstr));
            data = textscan(fid,'%f%f','Headerlines',1);
            fld = cell2mat(data);
            rdata     = fld(:,1);
            rdfdata   = fld(:,2);
            ax(i) = plot(rdata,rdfdata,'Color',pclr{i}, 'LineWidth', 2, 'MarkerSize', 8);
            
        end
        legendinfo{1} = 'Block-Block';
        legendinfo{2} = 'Block-Alter';
        legendinfo{3} = 'Alter-Block';
        legendinfo{4} = 'Alter-Alter';
        
        legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
        legend boxoff
        saveas(h1,sprintf('n%d_posions',nval),'png');
        
    end
end
%% Plot delU and delN

if delUplot == 1
    nvalsarr = [32,48,64,72];
    nmonfree = 30; nmongraft = 30; ngraft = 64;nsalt=510;nbackmons=10;
    ntotarr = nvalsarr*nmonfree + ngraft*(nmongraft+nbackmons) + 2*nsalt + ...
        abs(0.5*(ngraft-nvalsarr)*nmongraft);
    
    allUvals = importdata('./biascalc/delUint.txt');
    allNvals = importdata('./biascalc/delN.txt');
    
    varUext  = zeros(length(nvalsarr),4);
    meanUext = zeros(length(nvalsarr),4);
    
    varN  = zeros(length(nvalsarr),4);
    meanN = zeros(length(nvalsarr),4);
    
    
    %Cols 3-6 (U value std), 10-14 (U value mean)
    
    kvar = 3; kmean = 10;%See above comment
    for i = 1:4
        varUext(:,i) = allUvals.data(:,kvar);
        meanUext(:,i) = allUvals.data(:,kmean);
        kvar = kvar + 1;
        kmean = kmean + 1;
    end
    
    
    %Cols 3-6 (N value std), 10-14 (N value mean)
    
    kvar = 3; kmean = 10;%See above comment
    for i = 1:4
        varN(:,i) = allNvals.data(:,kvar);
        meanN(:,i) = allNvals.data(:,kmean);
        kvar = kvar + 1;
        kmean = kmean + 1;
    end
    
    varUext = varUext.*ntotarr;
    
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$N_{pa}/N_{pc}$','FontSize',20,'Interpreter','Latex')
    ylabel('$\langle \Delta U^{ext} \rangle$ ($k_B T$)','FontSize',20,'Interpreter','Latex')
    for i = 1:4
        errorbar(nvalsarr/ngraft,meanUext(:,i),varUext(:,i),'color',pclr{i},'LineWidth',2, ...
            'LineStyle',lsty{3},'Marker',msty{i},'MarkerSize',8,'MarkerFaceColor',pclr{i})
    end
    
    legendinfo{1} = 'Block-Block';
    legendinfo{2} = 'Block-Alter';
    legendinfo{3} = 'Alter-Block';
    legendinfo{4} = 'Alter-Alter';
    
    legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
    legend boxoff
    saveas(h1,'deltaUext_witherr','png');
    
    
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$N_{pa}/N_{pc}$','FontSize',20,'Interpreter','Latex')
    ylabel('$\langle \Delta N \rangle$','FontSize',20,'Interpreter','Latex')
    for i = 1:4
        errorbar(nvalsarr/ngraft,meanN(:,i),varN(:,i),'color',pclr{i},'LineWidth',2, ...
            'LineStyle',lsty{3},'Marker',msty{i},'MarkerSize',8,'MarkerFaceColor',pclr{i})
    end
    
    legendinfo{1} = 'Block-Block';
    legendinfo{2} = 'Block-Alter';
    legendinfo{3} = 'Alter-Block';
    legendinfo{4} = 'Alter-Alter';
    
    legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
    legend boxoff
    saveas(h1,'deltaUext_witherr','png');
    
end

%% Plot free energy Profiles
% For a particular N

if delFzplot   
    
    nvalsarr = [32];
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$z$','FontSize',20,'Interpreter','Latex')
    ylabel('$\Delta F(z)$ ($k_B T$)','FontSize',20,'Interpreter','Latex')
    for nvals = 1:length(nvalsarr)
        
        for i = 1:4
            
            if i == 1
                dirstr = 'block_block';
            elseif i == 2
                dirstr = 'block_alter';
            elseif i == 3
                dirstr = 'alter_block';
            else
                dirstr = 'alter_alter';
            end
            
            fylename = sprintf('./whamout_all/whamnew/n_%d_%s/whamout.txt',nvalsarr(nvals),dirstr);
            fid = fopen(fylename,'r');
            free_energy = zeros(10,3);
            header = fgetl(fid);
            
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
            
            errorbar(free_energy(:,1),free_energy(:,2),free_energy(:,3),'color',pclr{i},'LineWidth',2, ...
                'LineStyle',lsty{2},'Marker',msty{i},'MarkerSize',8,'MarkerFaceColor',pclr{i})
            
        end
        
        fclose(fid);
        
    end
    
    legendinfo{1} = 'Block-Block';
    legendinfo{2} = 'Block-Alter';
    legendinfo{3} = 'Alter-Block';
    legendinfo{4} = 'Alter-Alter';
    xlim([3 60])
    
    legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','NorthWest')
    legend boxoff
    saveas(h1,'delta_Fz','png');
end

%% Histogram Plots

if histplot == 1
    nvalsarr = [32]; nvals = 1;
    winarr  = {'5.0';'8.0';'10.0';'12.0';'14.0';'18.0';'22.0';'26.0';'30.0';'32.0';'34.0';'38.0';'42.0';'46.0'};
    
    for i = 1:1
        
        h1 = figure;
        hold on
        box on
        set(gca,'FontSize',16)
        xlabel('$z$','FontSize',20,'Interpreter','Latex')
        ylabel('Probability','FontSize',20,'Interpreter','Latex')
        
        for wval = 1:length(winarr)
            
            if i == 1
                dirstr = 'block_block';
            elseif i == 2
                dirstr = 'block_alter';
            elseif i == 3
                dirstr = 'alter_block';
            else
                dirstr = 'alter_alter';
            end
            
            fylename = sprintf('./whamout_all/USbackup_colvar/n_%d/out.colvars.traj_%s_%s',...
                nvalsarr(nvals),dirstr,winarr{wval});
            fid = fopen(fylename,'r');
            trval = zeros(1000,2);
            tline = fgetl(fid);
            k = 1;
            
            while ~feof(fid)
                
                tline = fgetl(fid);
                strarr = strsplit(tline);
                
                if strcmp(strarr{1},'#')
                    continue;
                else
                    trval(k,1) = str2double(strarr{2});
                    trval(k,2) = str2double(strarr{3});
                    k = k + 1;
                end
                
            end
            
            histogram(trval(:,2),80)
            
        end
        
        fclose(fid);
        saveas(h1,'hist','png');
        
    end
end

%% Fig Charge density, q(z)
% Ref:densananew.m



nfree = 32;rbin = 1;nbins=lz/rbin;
ncntr = abs((nfree-ngraft)*nmonfree/2);
ntotal = (nfreearr + ngraft)*nmonfree + 2*nsalt + ncntr + nbase;
chargearr = [0;1;0;-1;1;-1];

nval = nfree;
if nval > ngraft
    cntrcharg = 1;
elseif nval == ngraft
    cntrcharg = 0;
elseif nval < ngraft
    cntrcharg = -1;
end
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
    fid = fopen(sprintf('./results_dens/results_%d_%s/dens.lammpstrj',nval,dirstr));
    data = textscan(fid,'%f%f%f%f%f%f%f%f','Headerlines',1);
    fld = cell2mat(data);
    rdata     = fld(:,1);
    pneutral_brush = fld(:,2)*ngraft*nmongraft/2*lz*area/nbins;
    pnegativ_brush = fld(:,3)*ngraft*nmongraft/2*lz*area/nbins;
    pneutral_free  = fld(:,4)*nfree*nmonfree/2*lz*area/nbins;
    ppositive_free = fld(:,5)*nfree*nmonfree/2*lz*area/nbins;
    pposions  = fld(:,6)*nsalt*lz*area/nbins;
    pnegions  = fld(:,7)*nsalt*lz*area/nbins;
    pcntrions = fld(:,8)*ncntr*lz*area/nbins;
    
    if nval ~= ngraft
        netcharge = chargearr(1)*pneutral_brush + chargearr(2)*pnegativ_brush + ...
            chargearr(3)*pneutral_free + chargearr(4)*ppositive_free + ...
            chargearr(5)*pposions + chargearr(6)*pnegions + cntrcharg*pcntrions;
    else
        netcharge = chargearr(1)*pneutral_brush + chargearr(2)*pnegativ_brush + ...
            chargearr(3)*pneutral_free + chargearr(4)*ppositive_free + ...
            chargearr(5)*pposions + chargearr(6)*pnegions;
    end
    
    plot(rdata/lz,netcharge, 'color', pclr{i}, 'LineWidth', 2)
    
end
legendinfo{1} = 'Block-Block';
legendinfo{2} = 'Block-Alter';
legendinfo{3} = 'Alter-Block';
legendinfo{4} = 'Alter-Alter';

legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
legend boxoff
saveas(h1,sprintf('n%d_integral_netcharge',nval),'png');


%% Plot Ion density Profiles

if iondens
    nfree = 32;
    for j = 1:length(nfree)
        nval = nfree(j);
        h1 = figure;
        hold on
        box on
        set(gca,'FontSize',16)
        xlabel('$r/L_z$','FontSize',20,'Interpreter','Latex')
        ylabel('$\rho^{-}(r)$','FontSize',20,'Interpreter','Latex')
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
            fid = fopen(sprintf('./n_%d/results_%d_%s/dens.lammpstrj',nval,nval,dirstr));
            data = textscan(fid,'%f%f%f%f%f%f%f%f','Headerlines',1);
            fld = cell2mat(data);
            rdata     = fld(:,1);
            pposions  = fld(:,6);
            pnegions  = fld(:,7);
            pcntrions = fld(:,8);
            ax(i) = plot(rdata/lz,pnegions, pclr{i}, 'LineWidth', 2, 'MarkerSize', 8);
            %         plot(rdata/lz,pnegions, nclr{i}, 'LineWidth', 2, 'MarkerSize', 8)
            %         plot(rdata/lz,pcntrions, cclr{i}, 'LineWidth', 2, 'MarkerSize', 8)
        end
        legendinfo{1} = 'Block-Block';
        legendinfo{2} = 'Block-Alter';
        legendinfo{3} = 'Alter-Block';
        legendinfo{4} = 'Alter-Alter';
        
        legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
        legend boxoff
        saveas(h1,sprintf('n%d_posions',nval),'png');
        
    end
    
    for j = 1:length(nfree)
        nval = nfree(j);
        h1 = figure;
        hold on
        box on
        set(gca,'FontSize',16)
        xlabel('$r/L_z$','FontSize',20,'Interpreter','Latex')
        ylabel('$\rho^{+}(r)$','FontSize',20,'Interpreter','Latex')
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
            fid = fopen(sprintf('./n_%d/results_%d_%s/dens.lammpstrj',nval,nval,dirstr));
            data = textscan(fid,'%f%f%f%f%f%f%f%f','Headerlines',1);
            fld = cell2mat(data);
            rdata     = fld(:,1);
            pposions  = fld(:,6);
            pnegions  = fld(:,7);
            pcntrions = fld(:,8);
            ax(i) = plot(rdata/lz,pposions, pclr{i}, 'LineWidth', 2, 'MarkerSize', 8);
            %         plot(rdata/lz,pnegions, nclr{i}, 'LineWidth', 2, 'MarkerSize', 8)
            %         plot(rdata/lz,pcntrions, cclr{i}, 'LineWidth', 2, 'MarkerSize', 8)
        end
        legendinfo{1} = 'Block-Block';
        legendinfo{2} = 'Block-Alter';
        legendinfo{3} = 'Alter-Block';
        legendinfo{4} = 'Alter-Alter';
        
        legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
        legend boxoff
        saveas(h1,sprintf('n%d_posions',nval),'png');
        
    end
    
end