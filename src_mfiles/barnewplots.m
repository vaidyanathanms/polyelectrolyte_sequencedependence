%% To analyze the density profiles

clear;
clc;
close all;
format long;


%% Plotting Bar Plots of Adsorbed CHAINS

c = ({'Bl-Bl';'Bl-Al';'Al-Bl';'Al-Al'});

h11 = figure;
hold on
box on
set(gca,'FontSize',16)
ylabel('$f_{ads}$ (Chain Method)','FontSize',20,'Interpreter','Latex')
barvals = zeros(4);
xval = 1:1:4;
for i = 1:4
    if i == 1
        strval = 'blbl';
    elseif i == 2
        strval = 'blal';
    elseif i == 3
        strval = 'albl';
    elseif i == 4
        strval = 'alal';
    end
    data = importdata(sprintf('./newresults/chain_%s.txt',strval));
    barvals(i) = mean(data(:,2));
end
bar(barvals);
newlabel = {'Bl-Bl','Bl-Al','Al-Bl','Al-Al'};
set(gca,'XTick',1:4,'XtickLabel',{'Bl-Bl','Bl-Al','Al-Bl','Al-Al'});
saveas(h11,'chainmethod','png')

%% By number of Adsorbed MONOMERS
n = 100;
h11 = figure;
hold on
box on
set(gca,'FontSize',16)
ylabel('$f_{ads}$ (Monomers)','FontSize',20,'Interpreter','Latex')
barvals = zeros(4);

for i = 1:4
    if i == 1
        strval = 'blbl';
    elseif i == 2
        strval = 'blal';
    elseif i == 3
        strval = 'albl';
    elseif i == 4
        strval = 'alal';
    end
    data = importdata(sprintf('./newresults/free_%s.txt',strval));
    barvals(i) = mean(data(:,2));
end
bar(barvals);
newlabel = {'Bl-Bl','Bl-Al','Al-Bl','Al-Al'};
set(gca,'XTick',1:4,'XtickLabel',{'Bl-Bl','Bl-Al','Al-Bl','Al-Al'});
saveas(h11,'monomermethod','png')

%% By number of Adsorbed MONOMERS
n = 100;
h11 = figure;
hold on
box on
set(gca,'FontSize',16)
ylabel('$f_{ads}$ ($\int \rho_f(r) dr$ )','FontSize',20,'Interpreter','Latex')
barvals = zeros(4);

data = importdata(sprintf('./integ_htcutoff_%d.dat',n));
bar(data.data(:,2));

newlabel = {'Bl-Bl','Bl-Al','Al-Bl','Al-Al'};
set(gca,'XTick',1:4,'XtickLabel',{'Bl-Bl','Bl-Al','Al-Bl','Al-Al'});
saveas(h11,'integralmethod','png')