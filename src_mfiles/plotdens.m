%% To analyze the density profiles

clear;
clc;
close all;
format long;


%% Color Scheme

green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17]; brown=[0.6 0.2 0];
lsty = {'-','--',':','-.'};
msty = {'d','s','o','x'};
nclr = {'--r','--b','--k','--m','--g'};
cclr = {'r*','b*','k*','m*','g*'};
pclr = {'m',brown,green,'k','b', gold};

%% Plot Graft Profiles


nfree = [72]; lz = 120;
ngrafts = 64;area = 53^2;rbin = 1;nbins=lz/rbin;
nbase = 10*ngrafts;
nmons = 30;
nsalt = 510;
ncharg_chain = 15;
ntotal = (nfree + ngrafts)*nmons + 2*nsalt + nbase + ncharg_chain*(nfree+ngrafts);


%% Plot data

for j = 1:length(nfree)
    
    nval = nfree(j);
    
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
            legendinfo{2} = 'Neutral Block PC';
            legendinfo{3} = 'Charged Block PC ';
            legendinfo{4} = 'Neutral Block PA';
            legendinfo{5} = 'Charged Block PA';
        
        else
            
            plot(rdata/lz,pconnect, 'Color', pclr{1}, 'LineWidth', 2, 'MarkerSize', 8, 'LineStyle',lsty{1});
            plot(rdata/lz,pneutcat, 'Color', pclr{2}, 'LineWidth', 2, 'MarkerSize', 8, 'LineStyle',lsty{1});
            plot(rdata/lz,pcharcat, 'Color', pclr{3}, 'LineWidth', 2, 'MarkerSize', 8, 'LineStyle',lsty{2});
            %plot(rdata/lz,pneutani, 'Color', pclr{4}, 'LineWidth', 2, 'MarkerSize', 8, 'LineStyle',lsty{3});
            plot(rdata/lz,pcharani, 'Color', pclr{5}, 'LineWidth', 2, 'MarkerSize', 8, 'LineStyle',lsty{1});
            
            legendinfo{1} = 'Spacer';
            legendinfo{2} = 'Neutral Block PC';
            legendinfo{3} = 'Charged Block PC ';
            %legendinfo{4} = 'Neutral Block PA';
            legendinfo{4} = 'Charged Block PA';
        
        end
        
        legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
        legend boxoff
        xlim([0 0.30])
        saveas(h1,sprintf('./../../all_figures/n_%d_%s',nval,dirstr),'png');
        clear legendinfo
        
    end
    
end


