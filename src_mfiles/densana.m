%% To analyze the density profiles

clear;
clc;
close all;
format long;

%% Plot Ion Density Profiles

nmons = [32;64;80]; lz = 120;
pclr = {'r','b','k','m','g'};
nclr = {'--r','--b','--k','--m','--g'};
cclr = {'r*','b*','k*','m*','g*'};
    
% Positive Ions
h1 = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('$r/L_z$','FontSize',20,'Interpreter','Latex')
ylabel('$\rho^{+}(r)$','FontSize',20,'Interpreter','Latex')

for i = 1:length(nmons)
    fid = fopen(sprintf('./densprof/bl_bl/dens_%d.txt',nmons(i)));
    data = textscan(fid,'%f%f%f%f','Headerlines',1);
    fld = cell2mat(data);
    rdata   = fld(:,1);
    possalt = fld(:,2);
    negsalt = fld(:,3);
    counter = fld(:,4);
    plot(rdata/lz, possalt, pclr{i}, 'LineWidth', 2);
    %plot(rdata/lz, negsalt, nclr{i}, 'LineWidth', 2)
    %plot(rdata/lz, counter, cclr{i}, 'LineWidth', 2, 'MarkerSize', 8)
    legendinfo{i} = ['$N_g/N_f$: ' num2str(nmons(i)/64)];
end

legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
legend boxoff
saveas(h1,'posions','png');

h2 = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('$r/L_z$','FontSize',20,'Interpreter','Latex')
ylabel('$\rho^{-}(r)$','FontSize',20,'Interpreter','Latex')

for i = 1:length(nmons)
    fid = fopen(sprintf('./densprof/bl_bl/dens_%d.txt',nmons(i)));
    data = textscan(fid,'%f%f%f%f','Headerlines',1);
    fld = cell2mat(data);
    rdata   = fld(:,1);
    possalt = fld(:,2);
    negsalt = fld(:,3);
    counter = fld(:,4);
    %plot(rdata/lz, possalt, pclr{i}, 'LineWidth', 2);
    plot(rdata/lz, negsalt, nclr{i}, 'LineWidth', 2)
    %plot(rdata/lz, counter, cclr{i}, 'LineWidth', 2, 'MarkerSize', 8)
    legendinfo{i} = ['$N_g/N_f$: ' num2str(nmons(i)/64)];
end

legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
legend boxoff
saveas(h2,'negions','png');

h3 = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('$r/L_z$','FontSize',20,'Interpreter','Latex')
ylabel('$\rho^{c}(r)$','FontSize',20,'Interpreter','Latex')

for i = 1:length(nmons)
    fid = fopen(sprintf('./densprof/bl_bl/dens_%d.txt',nmons(i)));
    data = textscan(fid,'%f%f%f%f','Headerlines',1);
    fld = cell2mat(data);
    rdata   = fld(:,1);
    possalt = fld(:,2);
    negsalt = fld(:,3);
    counter = fld(:,4);
    %plot(rdata/lz, possalt, pclr{i}, 'LineWidth', 2);
    %plot(rdata/lz, negsalt, nclr{i}, 'LineWidth', 2)
    plot(rdata/lz, counter, cclr{i}, 'LineWidth', 2, 'MarkerSize', 8)
    legendinfo{i} = ['$N_g/N_f$: ' num2str(nmons(i)/64)];
end

legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
legend boxoff
saveas(h3,'cntions','png');

%% Plot Polymer Density Profiles

% Graft

h4 = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('$r/L_z$','FontSize',20,'Interpreter','Latex')
ylabel('$\rho^{g}(r)$','FontSize',20,'Interpreter','Latex')

for i = 1:length(nmons)
    fid = fopen(sprintf('./densprof/bl_bl/grp_%d.txt',nmons(i)));
    data = textscan(fid,'%f%f%f','Headerlines',1);
    fld = cell2mat(data);
    rdata   = fld(:,1);
    pegraft = fld(:,2);
    pefree  = fld(:,3);
    plot(rdata/lz,pegraft, pclr{i}, 'LineWidth', 2, 'MarkerSize', 8)
    legendinfo{i} = ['$N_g/N_f$: ' num2str(nmons(i)/64)];
end

legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
legend boxoff
saveas(h4,'graft','png');

% Free

h5 = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('$r/L_z$','FontSize',20,'Interpreter','Latex')
ylabel('$\rho^{f}(r)$','FontSize',20,'Interpreter','Latex')

for i = 1:length(nmons)
    fid = fopen(sprintf('./densprof/bl_bl/grp_%d.txt',nmons(i)));
    data = textscan(fid,'%f%f%f','Headerlines',1);
    fld = cell2mat(data);
    rdata   = fld(:,1);
    pegraft = fld(:,2);
    pefree  = fld(:,3);
    plot(rdata/lz,pefree, nclr{i}, 'LineWidth', 2, 'MarkerSize', 8)
    legendinfo{i} = ['$N_g/N_f$: ' num2str(nmons(i)/64)];
end

legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
legend boxoff
saveas(h5,'free','png');

clear legendinfo
% Plotting different architecture

n = 100;
h11 = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('$r/L_z$','FontSize',20,'Interpreter','Latex')
ylabel('$\rho^{f}(r)$','FontSize',20,'Interpreter','Latex')

for i = 1:4
    fid = fopen(sprintf('./densprof/n%d/grp_%d.txt',n,i));
    data = textscan(fid,'%f%f%f','Headerlines',1);
    fld = cell2mat(data);
    rdata   = fld(:,1);
    pegraft = fld(:,2);
    pefree  = fld(:,3);
    plot(rdata/lz,pefree, nclr{i}, 'LineWidth', 2, 'MarkerSize', 8)
end

legendinfo{1} = 'Block-Block';
legendinfo{2} = 'Block-Alter';
legendinfo{3} = 'Alter-Block';
legendinfo{4} = 'Alter-Alter';

legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
legend boxoff
saveas(h11,sprintf('n%d_free',n),'png');


h12 = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('$r/L_z$','FontSize',20,'Interpreter','Latex')
ylabel('$\rho^{f}(r)$','FontSize',20,'Interpreter','Latex')

for i = 1:4
    fid = fopen(sprintf('./densprof/n%d/grp_%d.txt',n,i));
    data = textscan(fid,'%f%f%f','Headerlines',1);
    fld = cell2mat(data);
    rdata   = fld(:,1);
    pegraft = fld(:,2);
    pefree  = fld(:,3);
    plot(rdata/lz,pegraft, nclr{i}, 'LineWidth', 2, 'MarkerSize', 8)
end

legendinfo{1} = 'Block-Block';
legendinfo{2} = 'Block-Alter';
legendinfo{3} = 'Alter-Block';
legendinfo{4} = 'Alter-Alter';

legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
legend boxoff
saveas(h12,sprintf('n%d_graft',n),'png');

