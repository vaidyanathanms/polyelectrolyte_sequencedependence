%% To plot all the energy quantities with errorbar as just VARIANCE

%Extra NOTE for computing bias energies:
%Technically internal energy is sum of kinetic and potential and hence both
%should be utilized. However, in a NVT ensemble, kinetic should cancel out 
%at %constant temeperature. But I am assuming it does not due to 
%fluctuations. Adding the kinetic energy term should NOT do anything to the
%final result.

clear
clc
close all
format long


%% Input Data

maxpointsperblock = 500;

nbiasmols = 5;
nmonfree = 30; nmongraft = 30; ngraft = 64; nbackmons = 10;
nfreearr = [32;48;64;72;80];
nsalt = 510;
ncounter = nfreearr*nmonfree/2 + ngraft*nmongraft/2;
cutoff = '1.50'; lz = 120; area=53^2;
rhofree = nfreearr*30/(lz*area);

green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17];brown = [0.2 0 0];
pclr = {'m',brown,green,'k','m', gold};
lsty = {'-','--',':'};
msty = {'d','s','o','x'};

freeener = 1;
biasU = 1;
unbiasU = 1;
adsfrac = 1;


%% Plot Free Energy

if freeener
    
    fylename = './../../all_txtfiles/expdeltaF_all.dat';
    if exist(fylename, 'file') == 2
        
        delF = zeros(length(nfreearr),4);
        k = 1;
        ffree = fopen(fylename,'r');
        fgetl(ffree);
        free_energy = zeros(length(nfreearr)*4,1);
        
        while ~feof(ffree) && k <= length(nfreearr)*4
            tline = fgetl(ffree);
            strarr = strsplit(tline);
            free_energy(k,1) = str2double(strarr{4});
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
        xlabel('$N_{pa}/N_{pc}$','FontSize',20,'Interpreter','Latex')
        ylabel('$\langle \Delta F \rangle$ ($k_B T$)','FontSize',20,'Interpreter','Latex')
        
        for i = 1:4
            plot(nfreearr/ngraft,delF(:,i),'color',pclr{i},'LineWidth',2,'LineStyle',...
                lsty{3},'Marker',msty{i},'MarkerSize',8,'MarkerFaceColor',pclr{i})
        end
        
        legendinfo{1} = 'Block-Block';
        legendinfo{2} = 'Block-Alter';
        legendinfo{3} = 'Alter-Block';
        legendinfo{4} = 'Alter-Alter';
        
        legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
        legend boxoff
        saveas(h1,'./../../all_figures/delFpow','png');
        
    end
end
%% Compute Bias Energy

if biasU
    avetemparr = zeros(length(nfreearr),4);
    aveenerarr = zeros(length(nfreearr),4);
    stdtemparr = zeros(length(nfreearr),4);
    stdenerarr = zeros(length(nfreearr),4);
    
    equilplusres = 1; % 1 = Both equilibrium and restart cycles are present in same file.
    neglectinit = 100; % neglect first so many cycles for equilibration
    
    if equilplusres == 1
        equilflag = 0; %Equilibrium end needs to be found
    else
        equilflag = 1; %Already in production stage
    end
    
    fout = fopen(sprintf('./../../all_txtfiles/AvgBiasEnergy_nbias_%d.dat',nbiasmols),'w');
    fprintf(fout,'%s\t%s\t%s\t%s\t%s\t%s\t%s\n','N_{pa}','Arch','N_tot','BiasedAveTemp(Intensive)',...
        'BiasedAveEner(Intensive)','BiasedAveTemp(Extensive)','BiasedAveEner(Extensive)');
    
    for ncnt = 1:length(nfreearr)
        
        nval = nfreearr(ncnt);
        nparticles = nval*nmonfree + ngraft*(nmongraft+nbackmons) + ...
            nsalt*2 + ncounter(ncnt);
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
            
            fylename = sprintf('./../../biascalc/ouputBiasEnergy_n_%d_%s.dat',...
                nval,dirstr);
            fin = fopen(fylename);
            
            fprintf('Analyzing Bias cycles for %d\t%s\n',nval,dirstr)
            
            % If both equilibrium and production cycles are present,
            % read until equilibrium cycles are over.
            
            if equilplusres == 1
                
                while equilflag == 0 && ~feof(fin)
                    
                    tline = strtrim(fgetl(fin));
                    strarr = strsplit(tline);
                    
                    if strcmp(strarr{1},'Step') %Found beginning of equil cycle
                        
                        tline = strtrim(fgetl(fin));
                        strarr = strsplit(tline);
                        
                        while ~strcmp(strarr{1},'Step') %Read until end of equil cycle
                            
                            tline = strtrim(fgetl(fin));
                            strarr = strsplit(tline);
                            
                        end
                        
                        equilflag = 1;
                        
                    end
                    
                end
                
            end
            
            if equilplusres == 1 && equilflag ~= 1
                fprintf('File Ended before finding production cycles for \n');
                fprintf('n=%d%s',nval, dirstr)
                continue;
            end
            
            while ~feof(fin)
                
                if strcmp(strarr{1},'Step')
                    
                    str2arr = strsplit(tline);
                    lenstr = length(str2arr);
                    
                    encol = -1; tempcol = -1;
                    % Extract columns for total energy and temperature
                    for sval = 1:lenstr
                        if strcmp(str2arr{sval},'TotEng')
                            encol = sval;
                        elseif strcmp(str2arr{sval},'Temp')
                            tempcol = sval;
                        end
                    end
                    
                    if encol == -1 || tempcol == -1
                        fprintf('Temp/Energy Column Flags = %d\t%d for n=%d%s',...
                            tempcol, encol, nval, dirstr)
                        break;
                    end
                    
                    temparr = zeros(100,1);
                    enerarr = zeros(100,1);
                    
                    sumtemp = 0; sumener = 0; cntvals = 0;
                    
                    % loop over neglectinit lines for equilibration
                    
                    for neglline = 1:neglectinit
                        if feof(fin)
                            fprintf('File Ended before finding production cycles for \n');
                            fprintf('n=%d%s',nval, dirstr)
                        end
                        tline = fgetl(fin);
                    end
                    
                    % Perform loop until end of file
                    while ~feof(fin)
                        
                        tline = strtrim(fgetl(fin));
                        strarr = strsplit(tline);
                        
                        if all(ismember(strarr{1}, '0123456789+-.eEdD'))
                            cntvals = cntvals + 1;
                            temparr(cntvals,1) = str2double(strarr{tempcol});
                            enerarr(cntvals,1) = str2double(strarr{encol});
                            sumtemp = sumtemp + temparr(cntvals);
                            sumener = sumener + enerarr(cntvals);
                        end
                        
                    end
                    %Technically internal energy is sum of kinetic and
                    %potential and hence both should be utilized. However,
                    %in a NVT ensemble, kinetic should cancel out at
                    %constant temeperature. But I am assuming it does not
                    %due to fluctuations. Adding the kinetic energy term
                    %should NOT do anything to the final result.
                    avetemparr(ncnt,i) = sumtemp/cntvals;
                    aveenerarr(ncnt,i) = sumener/cntvals;% - 1.5*avetemparr(ncnt,i);
                    
                    [tempvar, Tnblocksize] = blockave(temparr);
                    [enervar, Enblocksize] = blockave(enerarr);
                    
                    
                    fprintf(fout,'%d\t%s\t%d\t %16.8f\t %16.8f\t %16.8f\t %16.8f\n',...
                        nval,dirstr,nparticles,avetemparr(ncnt,i),aveenerarr(ncnt,i),...
                        avetemparr(ncnt,i)*nparticles,aveenerarr(ncnt,i)*nparticles);
                    
                    fvar = fopen(sprintf('./../../err_data/variancefileNoBias_n_%d_%s.dat',nval,dirstr),'w');
                    fprintf(fvar,'%s\t%s\t%s\t%s\t%s\t%s\n','EnergySteps','EnergyVar',...
                        'sqrt(enervar/N)','TempSteps','TempVar','sqrt(tempvar/N)');
                    
                    vartemp = var(temparr); varener = var(enerarr);
                    lentemp = length(temparr);
                    lenener = length(enerarr);
                    for varcnt = 1:maxpointsperblock
                        Tnchunks = floor(lentemp/varcnt);
                        Tnormstd = sqrt(tempvar(varcnt))/sqrt(Tnchunks);
                        Enchunks = floor(lenener/varcnt);
                        Enormstd = sqrt(enervar(varcnt))/sqrt(Enchunks);
                        fprintf(fvar,'%g\t%g\t%g\t%g\t%g\t%g\n',Enblocksize(varcnt), ...
                            enervar(varcnt),Enormstd,Tnblocksize(varcnt),tempvar(varcnt),Tnormstd);
                    end
                    stdenerarr(ncnt,i) = Enormstd; stdtemparr(ncnt,i) = Tnormstd;
                    fclose(fvar);
                    
                else
                    
                    tline = strtrim(fgetl(fin));
                    strarr = strsplit(tline);
                    
                end
                
            end
            
            fclose(fin);
            equilflag = 0;
            
        end
        
    end
    
    fin_biasaveener = aveenerarr;
    fin_biasstd = stdenerarr;
    
    %plot U-bias
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$N_{pa}/N_{pc}$','FontSize',20,'Interpreter','Latex')
    ylabel('$\langle U_{bias}^{int} \rangle $','FontSize',20,'Interpreter','Latex')
    
    for i = 1:4
        errorbar(nfreearr/ngraft,fin_biasaveener(:,i),fin_biasstd(:,i),'color',pclr{i},'LineWidth',2, ...
            'LineStyle',lsty{3},'Marker',msty{i},'MarkerSize',8,'MarkerFaceColor',pclr{i})
    end
    
    legendinfo{1} = 'Block-Block';
    legendinfo{2} = 'Block-Alter';
    legendinfo{3} = 'Alter-Block';
    legendinfo{4} = 'Alter-Alter';
    
    legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
    legend boxoff
    saveas(h1,'./../../all_figures/UbiasIntensive','png');
end

%% Compute Unbias Energy

if unbiasU
    avetemparr = zeros(length(nfreearr),4);
    aveenerarr = zeros(length(nfreearr),4);
    stdtemparr = zeros(length(nfreearr),4);
    stdenerarr = zeros(length(nfreearr),4);
    
    equilplusres = 0; % 1 = Both equilibrium and restart cycles are present in same file.
    
    if equilplusres == 1
        equilflag = 0; %Equilibrium end needs to be found
        neglectinit = 1250;
    else
        equilflag = 0; %Already in production stage
        neglectinit = 0;
    end
    
    fout = fopen('./../../all_txtfiles/AvgUnBiasEnergy_nbias.dat','w');
    fprintf(fout,'%s\t%s\t%s\t%s\t%s\t%s\t%s\n','N_{pa}','Arch','N_tot','UnBiasedAveTemp(Intensive)',...
        'UnBiasedAveEner(Intensive)','UnBiasedAveTemp(Extensive)','UnBiasedAveEner(Extensive)');
    
    for ncnt = 1:length(nfreearr)
        
        nval = nfreearr(ncnt);
        nparticles = nval*nmonfree + ngraft*(nmongraft+nbackmons) + ...
            nsalt*2 + ncounter(ncnt);
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
            
            fylename = sprintf('./../../biascalc/ouputNoBiasEnergy_n_%d_%s.dat',...
                nval,dirstr);
            fin = fopen(fylename);
            
            fprintf('Analyzing Unbias cycles for %d\t%s\n',nval,dirstr)
            
            % If both equilibrium and production cycles are present,
            % read until equilibrium cycles are over.
            
            if equilplusres == 1
                
                while equilflag == 0 && ~feof(fin)
                    
                    tline = strtrim(fgetl(fin));
                    strarr = strsplit(tline);
                    
                    if strcmp(strarr{1},'Step') %Found beginning of equil cycle
                        
                        tline = strtrim(fgetl(fin));
                        strarr = strsplit(tline);
                        
                        while ~strcmp(strarr{1},'Step') %Read until end of equil cycle
                            
                            tline = strtrim(fgetl(fin));
                            strarr = strsplit(tline);
                            
                        end
                        
                        equilflag = 1;
                        
                    end
                    
                end
                
            end
            
            if equilplusres == 1 && equilflag ~= 1
                fprintf('File Ended before finding production cycles for \n');
                fprintf('n=%d%s',nval, dirstr)
                continue;
            end
            
            while ~feof(fin)
                
                if strcmp(strarr{1},'Step')
                    
                    str2arr = strsplit(tline);
                    lenstr = length(str2arr);
                    
                    encol = -1; tempcol = -1;
                    % Extract columns for total energy and temperature
                    for sval = 1:lenstr
                        if strcmp(str2arr{sval},'TotEng')
                            encol = sval;
                        elseif strcmp(str2arr{sval},'Temp')
                            tempcol = sval;
                        end
                    end
                    
                    if encol == -1 || tempcol == -1
                        fprintf('Temp/Energy Column Flags = %d\t%d for n=%d%s',...
                            tempcol, encol, nval, dirstr)
                        break;
                    end
                    
                    temparr = zeros(100,1);
                    enerarr = zeros(100,1);
                    
                    sumtemp = 0; sumener = 0; cntvals = 0;
                    
                    % loop over neglectinit lines for equilibration
                    
                    for neglline = 1:neglectinit
                        if feof(fin)
                            fprintf('File Ended before finding production cycles for \n');
                            fprintf('n=%d%s',nval, dirstr)
                        end
                        tline = fgetl(fin);
                    end
                    
                    % Perform loop until end of file
                    while ~feof(fin)
                        
                        tline = strtrim(fgetl(fin));
                        if isempty(tline)
                            continue;
                        end
                        strarr = strsplit(tline);
                        
                        if all(ismember(strarr{1}, '0123456789+-.eEdD'))
                            cntvals = cntvals + 1;
                            temparr(cntvals,1) = str2double(strarr{tempcol});
                            enerarr(cntvals,1) = str2double(strarr{encol});
                            sumtemp = sumtemp + temparr(cntvals);
                            sumener = sumener + enerarr(cntvals);
                        end
                        
                    end
                    %Technically internal energy is sum of kinetic and
                    %potential and hence both should be utilized. However,
                    %in a NVT ensemble, kinetic should cancel out at
                    %constant temeperature. But I am assuming it does not
                    %due to fluctuations. Adding the kinetic energy term
                    %should NOT do anything to the final result.
                    avetemparr(ncnt,i) = sumtemp/cntvals;
                    aveenerarr(ncnt,i) = sumener/cntvals;% - 1.5*avetemparr(ncnt,i);
                    
                    [tempvar, Tnblocksize] = blockave(temparr);
                    [enervar, Enblocksize] = blockave(enerarr);
                    
                    
                    fprintf(fout,'%d\t%s\t%d\t %16.8f\t %16.8f\t %16.8f\t %16.8f\n',...
                        nval,dirstr,nparticles,avetemparr(ncnt,i),aveenerarr(ncnt,i),...
                        avetemparr(ncnt,i)*nbiasmols*nmonfree,aveenerarr(ncnt,i)*nbiasmols*nmonfree);
                    
                    fvar = fopen(sprintf('./../../err_data/variancefileNoBias_n_%d_%s.dat',nval,dirstr),'w');
                    fprintf(fvar,'%s\t%s\t%s\t%s\t%s\t%s\n','EnergySteps','EnergyVar',...
                        'sqrt(enervar/N)','TempSteps','TempVar','sqrt(tempvar/N)');
                    
                    vartemp = var(temparr); varener = var(enerarr);
                    lentemp = length(temparr);
                    lenener = length(enerarr);
                    for varcnt = 1:maxpointsperblock
                        Tnchunks = floor(lentemp/varcnt);
                        Tnormstd = sqrt(tempvar(varcnt))/sqrt(Tnchunks);
                        Enchunks = floor(lenener/varcnt);
                        Enormstd = sqrt(enervar(varcnt))/sqrt(Enchunks);
                        fprintf(fvar,'%g\t%g\t%g\t%g\t%g\t%g\n',Enblocksize(varcnt), ...
                            enervar(varcnt),Enormstd,Tnblocksize(varcnt),tempvar(varcnt),Tnormstd);
                    end
                    stdenerarr(ncnt,i) = Enormstd; stdtemparr(ncnt,i) = Tnormstd;
                    fclose(fvar);
                    
                else
                    
                    tline = strtrim(fgetl(fin));
                    strarr = strsplit(tline);
                    
                end
                
            end
            
            fclose(fin);
            equilflag = 0;
            
        end
        
    end
    fclose(fout);
    fin_Unbiasaveener = aveenerarr;
    fin_Unbiasstd = stdenerarr;
    
    %plot U-Unbias
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$N_{pa}/N_{pc}$','FontSize',20,'Interpreter','Latex')
    ylabel('$\langle U_{Nobias}^{int} \rangle$','FontSize',20,'Interpreter','Latex')
    
    for i = 1:4
        errorbar(nfreearr/ngraft,fin_Unbiasaveener(:,i),fin_Unbiasstd(:,i),'color',pclr{i},'LineWidth',2, ...
            'LineStyle',lsty{3},'Marker',msty{i},'MarkerSize',8,'MarkerFaceColor',pclr{i})
    end
    
    legendinfo{1} = 'Block-Block';
    legendinfo{2} = 'Block-Alter';
    legendinfo{3} = 'Alter-Block';
    legendinfo{4} = 'Alter-Alter';
    
    legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
    legend boxoff
    saveas(h1,'./../../all_figures/UNobias_intensive','png');
end
%% Analyze Delta N from Biased and Unbiased Simulations

if adsfrac
    fout = fopen(sprintf('./../../all_txtfiles/DeltaN_%s.dat',cutoff),'w');
    fprintf(fout,'%s\t%s\t%s\t%s\t%s\t%s\n','N_{pa}','Arch','N_bias','N_unbias',...
        'varNbias','DeltaN');
    
    avgN_bias = zeros(length(nfreearr),4);
    varN_bias = zeros(length(nfreearr),4);
    avgN_No_bias = zeros(length(nfreearr),4);
    DeltaN = zeros(length(nfreearr),4);
    
    for ncnt = 1:length(nfreearr)
        
        nval = nfreearr(ncnt);
        nparticles = nval*nmonfree + ngraft*(nmongraft+nbackmons) + ...
            nsalt*2 + ncounter(ncnt);
        
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
            
            fylename = sprintf('./../../biascalc/bias_adsfracchain_%s_n_%d_%s.dat',...
                cutoff,nval,dirstr);
            biasdata = importdata(fylename);
            
            fylename = sprintf('./../../biascalc/Nobias_adsfracchain_%s_n_%d_%s.dat',...
                cutoff,nval,dirstr);
            No_biasdata = importdata(fylename);
            
            avgN_bias(ncnt,i) = mean(biasdata(:,2));
            avgN_No_bias(ncnt,i) = mean(No_biasdata(:,2));
            DeltaN(ncnt,i) = avgN_No_bias(ncnt,i) - avgN_bias(ncnt,i);
            
            varN_bias(ncnt,i) = var(biasdata(:,2));
            
            fprintf(fout,'%d\t%s\t%g\t%g\t%g\t%g\n',nval,dirstr,avgN_bias(ncnt,i), ...
                avgN_No_bias(ncnt,i),varN_bias(ncnt,i),DeltaN(ncnt,i));
            
        end
        
    end
    fclose(fout);
    fin_delN = DeltaN;
    
    
    %plot U-Unbias
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$N_{pa}/N_{pc}$','FontSize',20,'Interpreter','Latex')
    ylabel('$\Delta N$','FontSize',20,'Interpreter','Latex')
    
    for i = 1:4
        errorbar(nfreearr/ngraft,fin_delN(:,i),varN_bias(:,1),'color',pclr{i},'LineWidth',2, ...
            'LineStyle',lsty{3},'Marker',msty{i},'MarkerSize',8,'MarkerFaceColor',pclr{i})
    end
    
    legendinfo{1} = 'Block-Block';
    legendinfo{2} = 'Block-Alter';
    legendinfo{3} = 'Alter-Block';
    legendinfo{4} = 'Alter-Alter';
    
    legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
    legend boxoff
    saveas(h1,'./../../all_figures/deltaN','png');
end
%% Plot Delta U and DeltaU/DeltaN

% Create a dummy particle array
partarr = zeros(length(nfreearr),4);
for ncnt = 1:length(nfreearr)
    nval = nfreearr(ncnt);
    for i = 1:4
        partarr(ncnt,i) = nval*nmonfree + ngraft*(nmongraft+nbackmons) + ...
            nsalt*2 + ncounter(ncnt);
    end
end

h1 = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('$N_{pa}/N_{pc}$','FontSize',20,'Interpreter','Latex')
ylabel('$\langle \Delta U^{ext} \rangle$','FontSize',20,'Interpreter','Latex')

for i = 1:4
    nparticles = partarr(:,i);
    varU = sqrt(nparticles).*sqrt(fin_biasstd.^2 + fin_Unbiasstd.^2);
    errorbar(nfreearr/ngraft,nparticles.*(fin_biasaveener(:,i)-fin_Unbiasaveener(:,i)),varU(:,i),...
        'color',pclr{i},'LineWidth',2,'LineStyle',lsty{3},'Marker',msty{i},...
        'MarkerSize',8,'MarkerFaceColor',pclr{i})
end

legendinfo{1} = 'Block-Block';
legendinfo{2} = 'Block-Alter';
legendinfo{3} = 'Alter-Block';
legendinfo{4} = 'Alter-Alter';

legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
legend boxoff
saveas(h1,'./../../all_figures/deltaU_extensive','png');

h1 = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('$N_{pa}/N_{pc}$','FontSize',20,'Interpreter','Latex')
ylabel('$\langle \Delta U \rangle$ ($k_B T$)','FontSize',20,'Interpreter','Latex') %actually delU/delN 
%in paper it is Delta U

for i = 1:4
    nparticles = partarr(:,i);
    varU = sqrt(nparticles).*sqrt(fin_biasstd.^2 + fin_Unbiasstd.^2);
    delU = partarr(:,i).*(fin_Unbiasaveener(:,i)-fin_biasaveener(:,i))./fin_delN(:,i);
    errorbar(nfreearr/ngraft,delU,varU(:,i),'color',pclr{i},'LineWidth',2,...
        'LineStyle',lsty{3},'Marker',msty{i},'MarkerSize',8,'MarkerFaceColor',pclr{i})
end

legendinfo{1} = 'Block-Block';
legendinfo{2} = 'Block-Alter';
legendinfo{3} = 'Alter-Block';
legendinfo{4} = 'Alter-Alter';

legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
legend boxoff
saveas(h1,'./../../all_figures/deltaUoverDeltaN','png');

%% Plot Delta S from delta U and Delta F

%load deltaF file

fylename = './../../all_txtfiles/expdeltaF_all.dat';
fsout = fopen('./../../all_txtfiles/AllEner.dat','w');
fprintf(fsout,'Number of brush/backbone/free mons:\t%d\t%d\t%d\n',nmonfree,nmongraft,nbackmons);
fprintf(fsout,'Number of graft chains/salt:\t%d\t%d\n',ngraft,nsalt);
fprintf(fsout,'%s\t%s\t%s\t%s\t%s\t%s\t%s\n','Nfree/Ngraft','Arch','deltaF','deltaU(int)','deltaN','delUbydelN','delS');
if exist(fylename, 'file') == 2
    
    delF = zeros(length(nfreearr),4);
    errF = zeros(length(nfreearr),4);
    k = 1;
    ffree = fopen(fylename,'r');
    header = fgetl(ffree);
    free_energy = zeros(length(nfreearr)*4,2);
    
    while ~feof(ffree) && k <= length(nfreearr)*4
        tline = fgetl(ffree);
        strarr = strsplit(tline);
        free_energy(k,1) = str2double(strarr{4});
        free_energy(k,2) = str2double(strarr{5});
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
    delU = zeros(length(nfreearr),4);
    delS = zeros(length(nfreearr),4);
    for i = 1:4
        nparticles = partarr(:,i);
        errorfromF = errF(:,i);
        varU = sqrt(nparticles.*(fin_biasstd.^2 + fin_Unbiasstd.^2)+ errorfromF.^2);
        delU(:,i) = partarr(:,i).*(fin_Unbiasaveener(:,i)-fin_biasaveener(:,i))./fin_delN(:,i);
        delS(:,i) = delU(:,i) - delF(:,i);
        errorbar(nfreearr/ngraft,delS(:,i),varU(:,i),'color',pclr{i},'LineWidth',2,...
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
                dirstr,delF(j,i),fin_biasaveener(j,i),fin_delN(j,i),delU(j,i),delS(j,i));
        end
    end
    legendinfo{1} = 'Block-Block';
    legendinfo{2} = 'Block-Alter';
    legendinfo{3} = 'Alter-Block';
    legendinfo{4} = 'Alter-Alter';
    
    legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
    legend boxoff
    saveas(h1,'./../../all_figures/deltaS','png');
    
    
else
    fprintf('%s not found',fylename);
    error('Cannot compute Delta S');
end
fclose(fsout);

