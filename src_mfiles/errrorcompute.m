%% To compute error bars in all quantities for energy analysis

clear
clc
close all
format long

%% Input Data

maxpointsperblock = 1000;

nbiasmols = 5;
nmonfree = 30; nmongraft = 30; ngraft = 64; nbackmons = 10;
nfreearr = [72];
nsalt = 510;
ncounter = abs(nfreearr*nmonfree/2 - ngraft*nmongraft/2);
cutoff = '1.50'; lz = 120; area=53^2;
rhofree = nfreearr*30/(lz*area);

green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17];
pclr = {'r','b',green,'k','m', gold};
lsty = {'-','--',':'};
msty = {'d','s','o','x'};


%%  Compute error bar in N

for ncnt = 1:length(nfreearr)
    
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$N_f/N_{pc}$','FontSize',20,'Interpreter','Latex')
    ylabel('$N_{bias}$','FontSize',20,'Interpreter','Latex')
    
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
        
        fylename = sprintf('./biasall_results/bias_adsfracchain_%s_n_%d_%s.dat',...
            cutoff,nval,dirstr);
        biasdata = importdata(fylename);
        
        [bvar,svar] = blockave(biasdata(:,2));
        plot(svar,bvar.*sqrt(svar/length(biasdata(:,2))))
        
    end
    
    legendinfo{1} = 'Block-Block';
    legendinfo{2} = 'Alter-Block';
    legendinfo{3} = 'Block-Alter';
    legendinfo{4} = 'Alter-Alter';
    legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
    legend boxoff
    saveas(h1,sprintf('stddev_biasdelN_n%d.png',nval));
    
end

for ncnt = 1:length(nfreearr)
    
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$N_{pa}/N_{pc}$','FontSize',20,'Interpreter','Latex')
    ylabel('$N_{unbias}$','FontSize',20,'Interpreter','Latex')
    
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
        
        
        fylename = sprintf('./biasall_results/Nobias_adsfracchain_%s_n_%d_%s.dat',...
            cutoff,nval,dirstr);
        No_biasdata = importdata(fylename);
        
        [bvar,svar] = blockave(No_biasdata(:,2));
        plot(svar,bvar.*sqrt(svar/length(No_biasdata(:,2))))
        
    end
    
    legendinfo{1} = 'Block-Block';
    legendinfo{2} = 'Alter-Block';
    legendinfo{3} = 'Block-Alter';
    legendinfo{4} = 'Alter-Alter';
    legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
    legend boxoff
    saveas(h1,sprintf('stddev_no_biasdelN_n%d.png',nval));
    
end

%% Compute Bias Energy

equilplusres = 1; % 1 = Both equilibrium and restart cycles are present in same file.

if equilplusres == 1
    equilflag = 0; %Equilibrium end needs to be found
else
    equilflag = 1; %Already in production stage
end

fout = fopen(sprintf('AvgBiasEnergy_nbias_%d.dat',nbiasmols),'w');
fprintf(fout,'%s\t%s\t%s\t%s\t%s\t%s\t%s\n','N_{pa}','Arch','N_tot','BiasedAveTemp(Intensive)',...
    'BiasedAveEner(Intensive)','BiasedAveTemp(Extensive)','BiasedAveEner(Extensive)');

for ncnt = 1:length(nfreearr)
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$N_{pa}/N_{pc}$','FontSize',20,'Interpreter','Latex')
    ylabel('$\Delta U_{bias}$','FontSize',20,'Interpreter','Latex')
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
        
        fylename = sprintf('./biasall_results/ouputBiasEnergy_n_%d_%s.dat',...
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
                
                
                [enervar, Enblocksize] = blockave(enerarr);
                plot(Enblocksize,enervar.*sqrt(Enblocksize/length(enerarr(:,1))))
                
            else
                
                tline = strtrim(fgetl(fin));
                strarr = strsplit(tline);
                
            end
            
        end
        
        fclose(fin);
        equilflag = 0;
        
    end
    
    
    legendinfo{1} = 'Block-Block';
    legendinfo{2} = 'Block-Alter';
    legendinfo{3} = 'Alter-Block';
    legendinfo{4} = 'Alter-Alter';
    legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
    legend boxoff
    saveas(h1,sprintf('stddev_UbiasIntensive_%d',nval),'png');
    
end



%% Compute Unbias Energy

equilplusres = 0; % 1 = Both equilibrium and restart cycles are present in same file.

if equilplusres == 1
    equilflag = 0; %Equilibrium end needs to be found
else
    equilflag = 1; %Already in production stage
end

for ncnt = 1:length(nfreearr)
    h1 = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('$N_{pa}/N_{pc}$','FontSize',20,'Interpreter','Latex')
    ylabel('$\Delta U_{unbias}$','FontSize',20,'Interpreter','Latex')
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
        
        fylename = sprintf('./biasall_results/ouputNoBiasEnergy_n_%d_%s.dat',...
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
                
                
                [enervar, Enblocksize] = blockave(enerarr);
                plot(Enblocksize,enervar.*sqrt(Enblocksize/length(enerarr(:,1))))
                
            else
                
                tline = strtrim(fgetl(fin));
                strarr = strsplit(tline);
                
            end
            
        end
        
        fclose(fin);
        equilflag = 0;
        
    end
    legendinfo{1} = 'Block-Block';
    legendinfo{2} = 'Block-Alter';
    legendinfo{3} = 'Alter-Block';
    legendinfo{4} = 'Alter-Alter';
    legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
    legend boxoff
    saveas(h1,sprintf('stddev_No_UbiasIntensive_%d',nval),'png');
end
