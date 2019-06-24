%% To analyze the density profiles

clear;
clc;
close all;
format long;

%% Color Scheme

green = [0 0.5 0.0]; gold = [0.9 0.75 0];
orange = [0.91 0.41 0.17]; brown=[0.6 0.2 0];
pclr = {'m',brown,green,'k','b', gold};
lsty = {'-','--',':'};
msty = {'d','s','o','x'};

%% Input data

nfree = [32;48;64;72]; lz = 120; 
ngrafts = 64;area = 53^2;rbin = 1;nbins=lz/rbin;
nbase = 10;
nmons = 30;
nsalt = 510;
ncharg_chain = 15;
ntotal = (nfree + ngrafts)*nmons + 2*nsalt + nbase*ngrafts + ncharg_chain*(nfree+ngrafts);

%% Plot Polymer Graft Profiles

for j = 1:length(nfree)
    
    nval = nfree(j);
    
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$r/L_z$','FontSize',20,'Interpreter','Latex')
    ylabel('$\rho_{f}(r)$','FontSize',20,'Interpreter','Latex')
    
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
        
        sprintf('./../../results_dens/results_%d_%s/grpdens.lammpstrj',nval,dirstr);
        fid = fopen(sprintf('./../../results_dens/results_%d_%s/grpdens.lammpstrj',nval,dirstr));
        
        data = textscan(fid,'%f%f%f','Headerlines',1);
        fld = cell2mat(data);
        
        rdata   = fld(:,1);
        pegraft = fld(:,2);
        pefree  = fld(:,3);
        
        plot(rdata/lz, pegraft, 'Color', pclr{i}, 'LineWidth', 2, 'MarkerSize', 8)
        
    end
    
    legendinfo{1} = 'Block-Block';
    legendinfo{2} = 'Block-Alter';
    legendinfo{3} = 'Alter-Block';
    legendinfo{4} = 'Alter-Alter';
    
    legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
    legend boxoff
    saveas(h1,sprintf('./../../all_figures/n%d_graft',nval),'png');
    
end

%% Plot Polymer free Profiles

for j = 1:length(nfree)
    
    nval = nfree(j);
    
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$r/L_z$','FontSize',20,'Interpreter','Latex')
    ylabel('$\rho_{f}(r)$','FontSize',20,'Interpreter','Latex')
    
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
        
        sprintf('./../../results_dens/results_%d_%s/grpdens.lammpstrj',nval,dirstr);
        fid = fopen(sprintf('./../../results_dens/results_%d_%s/grpdens.lammpstrj',nval,dirstr));
        data = textscan(fid,'%f%f%f','Headerlines',1);
        fld = cell2mat(data);
    
        rdata   = fld(:,1);
        pegraft = fld(:,2);
        pefree  = fld(:,3);
    
        plot(rdata/lz,pefree, 'Color', pclr{i}, 'LineWidth', 2, 'MarkerSize', 8)
    
    end
    
    legendinfo{1} = 'Block-Block';
    legendinfo{2} = 'Block-Alter';
    legendinfo{3} = 'Alter-Block';
    legendinfo{4} = 'Alter-Alter';
    
    legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
    legend boxoff
    saveas(h1,sprintf('./../../all_figures/n%d_free',nval),'png');
    
end

%% Plot Ion density Profiles

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
        
        fid = fopen(sprintf('./../../results_dens/results_%d_%s/dens.lammpstrj',nval,dirstr));
        data = textscan(fid,'%f%f%f%f%f%f%f%f','Headerlines',1);
        fld = cell2mat(data);
        
        rdata     = fld(:,1);
        pposions  = fld(:,7);
        pnegions  = fld(:,8);
        
        ax(i) = plot(rdata/lz,pposions, 'Color', pclr{i}, 'LineWidth', 2, 'MarkerSize', 8);
        
    end
    
    legendinfo{1} = 'Block-Block';
    legendinfo{2} = 'Block-Alter';
    legendinfo{3} = 'Alter-Block';
    legendinfo{4} = 'Alter-Alter';
    
    legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
    legend boxoff
    saveas(h1,sprintf('./../../all_figures/n%d_posions',nval),'png');
    
end

%% Plot Overall Charge

chargearr = [0;0;1;0;-1;1;-1];

for j = 1:length(nfree)

    nval = nfree(j);

    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$r/L_z$','FontSize',20,'Interpreter','Latex')
    ylabel('$Q(r)$','FontSize',20,'Interpreter','Latex')

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
        
        fid = fopen(sprintf('./../../results_dens/results_%d_%s/dens.lammpstrj',nval,dirstr));
        data = textscan(fid,'%f%f%f%f%f%f%f%f','Headerlines',1);
        fld = cell2mat(data);
        
        nposions = nsalt + nval*ncharg_chain;
        nnegions = nsalt + ngrafts*ncharg_chain;
        
        rdata     = fld(:,1);
        pbackbon_brush = fld(:,2)*ngrafts*nbase*lz*area/nbins;
        pneutral_brush = fld(:,3)*ngrafts*nmons/2*lz*area/nbins;
        pnegativ_brush = fld(:,4)*ngrafts*nmons/2*lz*area/nbins;
        pneutral_free  = fld(:,5)*nfree(j)*nmons/2*lz*area/nbins;
        ppositive_free = fld(:,6)*nfree(j)*nmons/2*lz*area/nbins;
        pposions  = fld(:,7)*nposions*lz*area/nbins;
        pnegions  = fld(:,8)*nnegions*lz*area/nbins;
        
        netcharge = chargearr(1)*pbackbon_brush + chargearr(2)*pneutral_brush ... 
            + chargearr(3)*pnegativ_brush + chargearr(4)*pneutral_free + ... 
            chargearr(5)*ppositive_free + chargearr(6)*pposions + chargearr(7)*pnegions;
        
        plot(rdata/lz,netcharge, 'Color', pclr{i}, 'LineWidth', 2, 'MarkerSize', 8)
    
    end
    
    legendinfo{1} = 'Block-Block';
    legendinfo{2} = 'Block-Alter';
    legendinfo{3} = 'Alter-Block';
    legendinfo{4} = 'Alter-Alter';
    
    legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
    legend boxoff
    saveas(h1,sprintf('./../../all_figures/n%d_netcharge',nval),'png');
    
end


%% Plot Integral Charge

chargearr = [0;0;1;0;-1;1;-1];


for j = 1:length(nfree)
    
    nval = nfree(j);
    
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$r/L_z$','FontSize',20,'Interpreter','Latex')
    ylabel('$ Q(r) $','FontSize',20,'Interpreter','Latex')
    
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
        
        fid = fopen(sprintf('./../../results_dens/results_%d_%s/dens.lammpstrj',nval,dirstr));
        data = textscan(fid,'%f%f%f%f%f%f%f%f','Headerlines',1);
        fld = cell2mat(data);
        
        nposions = nsalt + nval*ncharg_chain;
        nnegions = nsalt + ngrafts*ncharg_chain;
        
        rdata     = fld(:,1);
        pbackbon_brush = fld(:,2)*ngrafts*nbase*lz*area/nbins;
        pneutral_brush = fld(:,3)*ngrafts*nmons/2*lz*area/nbins;
        pnegativ_brush = fld(:,4)*ngrafts*nmons/2*lz*area/nbins;
        pneutral_free  = fld(:,5)*nfree(j)*nmons/2*lz*area/nbins;
        ppositive_free = fld(:,6)*nfree(j)*nmons/2*lz*area/nbins;
        pposions  = fld(:,7)*nposions*lz*area/nbins;
        pnegions  = fld(:,8)*nnegions*lz*area/nbins;
        
        netcharge = chargearr(1)*pbackbon_brush + chargearr(2)*pneutral_brush ... 
            + chargearr(3)*pnegativ_brush + chargearr(4)*pneutral_free + ... 
            chargearr(5)*ppositive_free + chargearr(6)*pposions + chargearr(7)*pnegions;
               
        intnet = zeros(length(netcharge),1);
        intnet(1,1) = 0.5*netcharge(1);
        rnet(1,1) = 0.5*(rdata(1));
        
        for k = 2:length(netcharge)
            intnet(k,1) = intnet(k-1,1) + 0.5*(rdata(k)-rdata(k-1))*(netcharge(k)+netcharge(k-1));
            rnet(k,1) = 0.5*(rdata(k-1,1)+rdata(k,1));
        end
        
        plot(rnet/lz,intnet, 'Color', pclr{i}, 'LineWidth', 2, 'MarkerSize', 8)
    
    end
    
    legendinfo{1} = 'Block-Block';
    legendinfo{2} = 'Block-Alter';
    legendinfo{3} = 'Alter-Block';
    legendinfo{4} = 'Alter-Alter';
    
    legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
    legend boxoff
    saveas(h1,sprintf('./../../all_figures/n%d_integral_netcharge',nval),'png');
    
end

