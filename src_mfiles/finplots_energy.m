%% To plot all the energy quantities with actual errorbar
clear
clc
close all
format long


%% Input Data

maxpointsperblock = 1000;

nbiasmols = 5;
nmonfree = 30; nmongraft = 30; ngraft = 64; nbackmons = 10;
nfreearr = [32;48;64;72;80;100];
nsalt = 510;
ncounter = abs(nfreearr*nmonfree/2 - ngraft*nmongraft/2);
cutoff = '1.50'; lz = 120; area=53^2;
rhofree = nfreearr*30/(lz*area);

green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17];
pclr = {'r','b',green,'k','m', gold};
lsty = {'-','--',':'};
msty = {'d','s','o','x'};

%% delta N with errorbar
dataN = importdata('ndata.txt');

%Cols 3-6 (N value std), Cols 10-13 (N value mean)

varN  = zeros(length(nfreearr),4);
meanN = zeros(length(nfreearr),4);

kvar = 3; kmean = 10; %See above comment
for i = 1:4
    varN(:,i) = dataN.data(:,kvar);
    meanN(:,i) = dataN.data(:,kmean);
    kvar = kvar + 1; kmean = kmean+1;
end

h1 = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('$N_f/N_g$','FontSize',20,'Interpreter','Latex')
ylabel('$\langle \Delta N \rangle$','FontSize',20,'Interpreter','Latex')
for i = 1:4
    errorbar(nfreearr/ngraft,meanN(:,i),varN(:,i),'color',pclr{i},'LineWidth',2, ...
        'LineStyle',lsty{3},'Marker',msty{i},'MarkerSize',8,'MarkerFaceColor',pclr{i})
end

legendinfo{1} = 'Block-Block';
legendinfo{2} = 'Block-Alter';
legendinfo{3} = 'Alter-Block';
legendinfo{4} = 'Alter-Alter';

legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
legend boxoff
saveas(h1,'deltaN_witherr','png');


%% Extensive internal energy with errorbar

dataU = importdata('Udata.txt');

%Cols 3-6 (N value std), Cols 10-13 (N value mean)

varU  = zeros(length(nfreearr),4);
meanU = zeros(length(nfreearr),4);
ntot  = dataU.data(:,2);

kvar = 3; kmean = 10; %See above comment
for i = 1:4
    varU(:,i) = ntot(:,1).*dataU.data(:,kvar);
    meanU(:,i) = dataU.data(:,kmean);
    kvar = kvar + 1; kmean = kmean+1;
end

h1 = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('$N_f/N_g$','FontSize',20,'Interpreter','Latex')
ylabel('$\langle \Delta U^{ext} \rangle$','FontSize',20,'Interpreter','Latex')
for i = 1:4
    errorbar(nfreearr/ngraft,meanU(:,i),varU(:,i),'color',pclr{i},'LineWidth',2, ...
        'LineStyle',lsty{3},'Marker',msty{i},'MarkerSize',8,'MarkerFaceColor',pclr{i})
end

legendinfo{1} = 'Block-Block';
legendinfo{2} = 'Block-Alter';
legendinfo{3} = 'Alter-Block';
legendinfo{4} = 'Alter-Alter';

legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
legend boxoff
saveas(h1,'deltaUext_witherr','png');

%% Intensive (delU/delN) with errorbar

dataUbyN = importdata('deludelN.txt');

%Cols 3-6 (N value std), Cols 10-13 (N value mean)

varUbyN  = zeros(length(nfreearr),4);
meanUbyN = zeros(length(nfreearr),4);

kvar = 3; kmean = 10; %See above comment
for i = 1:4
    varUbyN(:,i) = dataUbyN.data(:,kvar);
    meanUbyN(:,i) = dataUbyN.data(:,kmean);
    kvar = kvar + 1; kmean = kmean+1;
end

h1 = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('$N_f/N_g$','FontSize',20,'Interpreter','Latex')
ylabel('$\langle \Delta U \rangle$','FontSize',20,'Interpreter','Latex')
for i = 1:4
    errorbar(nfreearr/ngraft,meanUbyN(:,i),varUbyN(:,i),'color',pclr{i},'LineWidth',2, ...
        'LineStyle',lsty{3},'Marker',msty{i},'MarkerSize',8,'MarkerFaceColor',pclr{i})
end

legendinfo{1} = 'Block-Block';
legendinfo{2} = 'Block-Alter';
legendinfo{3} = 'Alter-Block';
legendinfo{4} = 'Alter-Alter';

legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
legend boxoff
saveas(h1,'deltaN_witherr','png');

%% entropy with errorbar

fylename = './../deltaF_all.dat';
fsout = fopen('AllEner.dat','w');
fprintf(fsout,'Number of brush/backbone/free mons:\t%d\t%d\t%d\n',nmonfree,nmongraft,nbackmons);
fprintf(fsout,'Number of graft chains/salt:\t%d\t%d\n',ngraft,nsalt);
fprintf(fsout,'%s\t%s\t%s\t%s\t%s\t%s\n','Nfree/Ngraft','Arch','deltaF','deltaN','delUbydelN','delS');
if exist(fylename, 'file') == 2
    
    delF = zeros(length(nfreearr),4);
    k = 1;
    ffree = fopen(fylename,'r');
    header = fgetl(ffree);
    free_energy = zeros(length(nfreearr)*4,1);
    
    while ~feof(ffree) && k <= length(nfreearr)*4
        tline = fgetl(ffree);
        strarr = strsplit(tline);
        free_energy(k,1) = str2double(strarr{3});
        k = k + 1;
    end
    fclose(ffree);
    
    if k ~= length(nfreearr)*4+1
        error('Unequal number of free energy columns%d\t%d\n',k,length(nfreearr)*4+1);
    end
    
    k = 1;
    for i = 1:length(nfreearr)
        for j = 1:4
            delF(i,j) = free_energy(k,1);
            k = k + 1;
        end
    end
    
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$N_f/N_g$','FontSize',20,'Interpreter','Latex')
    ylabel('$\langle T \Delta S \rangle$','FontSize',20,'Interpreter','Latex')
    
    for i = 1:4        
        delS = meanUbyN(:,i) - delF(:,i);
        errorbar(nfreearr/ngraft,delS,varUbyN(:,1),'color',pclr{i},'LineWidth',2, ...
            'LineStyle',lsty{3},'Marker',msty{i},'MarkerSize',8,'MarkerFaceColor',pclr{i})
        
        
        for j = 1:length(nfreearr)
            if i == 1
                dirstr = 'bl_bl';
            elseif i == 2
                dirstr = 'bl_al';
            elseif i == 3
                dirstr = 'al_bl';
            else
                dirstr = 'al_al';
            end
            fprintf(fsout,'%g\t%s\t%g\t%g\t%g\t%g\t%g\n',nfreearr(j)/ngraft,...
                dirstr,delF(j,i),meanN(j,i),meanUbyN(j,i),delS(j,1));
        end
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
fclose(fsout);
