%% To compute average adsorbed fraction of chains
% 1) Use this first to find the errors in data and then use the
% finplots_forpaper.m to plot using errorbars.
% 2) Change path to filename (search for keyword filename) and other places
% wherever necessary.

clear;
clc;
close all;
format long;

%% Inputs

nmonfree = 30; nmongraft = 30; ngraft = 64;
nfreearr = [16;32;48;64;72;80;100];
cutoff = '2.00'; lz = 120; area=53^2;
rhofree = nfreearr*30/(lz*area);


green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17];brown = [0.2 0 0];
pclr = {'m',brown,green,'k','m', gold};
lsty = {'-','--',':'};
msty = {'d','s','o','x'};

nadschain = zeros(length(nfreearr),4);
stddevon = 0;
%% Plot adsorbed fraction as a function of number of ADSORBED CHAINS
fout = fopen(sprintf('./../../all_txtfiles/adsorbed_chain_ave_rcut_%s.dat',cutoff),'w');
fprintf(fout,'%s\t%s\t%s\n','N_f','Arch','fraction');
fprintf('%s\t%s\t%s\n','N_f','Arch','fraction');

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
        
        filename = sprintf('./../../results_adsfrac/all_cutoff/results_%d_%s/adsfracchain_rcut_%s.lammpstrj',...
            nval,dirstr,cutoff);
        data = importdata(filename);
        if min(size(data)) == 0
            fprintf('No data for %d\t%s', nval, dirstr)
            continue;
        end
        nadschain(ncnt,i) = mean(data(:,3));
        fprintf(fout,'%d\t%s\t%g\n',nval,dirstr,nadschain(ncnt,i));
        fprintf('%d\t%s\t%g\n',nval,dirstr,nadschain(ncnt,i));
    end
end
fclose(fout);

h1 = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('$N_{pa}/N_{pc}$','FontSize',20,'Interpreter','Latex')
ylabel('$f_{ads}$','FontSize',20,'Interpreter','Latex')

plot(nfreearr/ngraft,nadschain(:,1),'color',pclr{1},'LineWidth',2,'LineStyle',lsty{3},'Marker',msty{1},'MarkerSize',8,'MarkerFaceColor',pclr{1})
plot(nfreearr/ngraft,nadschain(:,3),'color',pclr{3},'LineWidth',2,'LineStyle',lsty{3},'Marker',msty{3},'MarkerSize',8,'MarkerFaceColor',pclr{3})
plot(nfreearr/ngraft,nadschain(:,2),'color',pclr{2},'LineWidth',2,'LineStyle',lsty{3},'Marker',msty{2},'MarkerSize',8,'MarkerFaceColor',pclr{2})
plot(nfreearr/ngraft,nadschain(:,4),'color',pclr{4},'LineWidth',2,'LineStyle',lsty{3},'Marker',msty{4},'MarkerSize',8,'MarkerFaceColor',pclr{4})

legendinfo{1} = 'Block-Block';
legendinfo{2} = 'Alter-Block';
legendinfo{3} = 'Block-Alter';
legendinfo{4} = 'Alter-Alter';


%overlay y = x line

x = 0:0.1:1.1; y = x;
plot(x,y,'LineWidth',2,'Color',orange,'LineStyle','--')


legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
legend boxoff
saveas(h1,sprintf('./../../all_figures/figs_repo/adsorbchain_rcut_%s.png',cutoff));

%compute standard dev

if stddevon
    for ncnt = 1:length(nfreearr)
        nval = nfreearr(ncnt);
        h1 = figure;
        hold on
        box on
        set(gca,'FontSize',16)
        ylabel('BSE','FontSize',20,'Interpreter','Latex')
        xlabel('Time','FontSize',20,'Interpreter','Latex')
        %title(['$n =$ ' num2str(nval)], 'FontSize',20,'Interpreter','Latex')
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
            filename = sprintf('./../../results_adsfrac/results_%d_%s/adsfracchain_rcut_%s.lammpstrj',...
                nval,dirstr,cutoff);
            data = importdata(filename);
            nadschain(ncnt,i) = mean(data(:,3));
            [bvar,svar] = blockave(data(:,3));
            plot(svar,bvar.*sqrt(svar/length(data(:,3))))
        end
        legendinfo{1} = 'Block-Block';
        legendinfo{2} = 'Alter-Block';
        legendinfo{3} = 'Block-Alter';
        legendinfo{4} = 'Alter-Alter';
        legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
        legend boxoff
        saveas(h1,sprintf('./../../err_data/figs_err/stddev_adsfrac_n%d.fig',nval));
        saveas(h1,sprintf('./../../err_data/figs_err/stddev_adsfrac_n%d.png',nval));
    end
end




%% Plot adsorbed fraction as a function of number of ADSORBED MONOMERS
fout = fopen(sprintf('./../../all_txtfiles/adsorbed_mon_ave_rcut_%s.dat',cutoff),'w');
fprintf(fout,'%s\t%s\t%s\n','N_f','Arch','fraction');

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
        
        filename = sprintf('./../../results_adsfrac/results_%d_%s/adsfracmon_rcut_%s.lammpstrj',...
            nval,dirstr,cutoff);
        data = importdata(filename);
        nadschain(ncnt,i) = nval*mean(data(:,3))/ngraft;
        fprintf(fout,'%d\t%s\t%g\n',nval,dirstr,nadschain(ncnt,i));
    end
end
fclose(fout);

h1 = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('$N_f/N_g$','FontSize',20,'Interpreter','Latex')
ylabel('$f$','FontSize',20,'Interpreter','Latex')

plot(nfreearr/ngraft,nadschain(:,1),'color',pclr{1},'LineWidth',2,'LineStyle',lsty{3},'Marker',msty{1},'MarkerSize',8,'MarkerFaceColor',pclr{1})
plot(nfreearr/ngraft,nadschain(:,3),'color',pclr{3},'LineWidth',2,'LineStyle',lsty{3},'Marker',msty{3},'MarkerSize',8,'MarkerFaceColor',pclr{3})
plot(nfreearr/ngraft,nadschain(:,2),'color',pclr{2},'LineWidth',2,'LineStyle',lsty{3},'Marker',msty{2},'MarkerSize',8,'MarkerFaceColor',pclr{2})
plot(nfreearr/ngraft,nadschain(:,4),'color',pclr{4},'LineWidth',2,'LineStyle',lsty{3},'Marker',msty{4},'MarkerSize',8,'MarkerFaceColor',pclr{4})

legendinfo{1} = 'Block-Block';
legendinfo{2} = 'Alter-Block';
legendinfo{3} = 'Block-Alter';
legendinfo{4} = 'Alter-Alter';

legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
legend boxoff
saveas(h1,sprintf('./../../all_figures/figs_repo/adsorbmon_rcut_%s.png',cutoff));