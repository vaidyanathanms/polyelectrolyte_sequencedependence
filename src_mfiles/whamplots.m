clc;
clear;
close all;
format long;


%% Input Flags

histflag = 0;
delFvsfofz = 0;

%% Input Data

nvalsarr = [32,48,64,72,80]; ngraft = 64;
nchargedmons = 15; nsalt = 510;
diff_energy = zeros(length(nvalsarr),4);
err_energy = zeros(length(nvalsarr),4);
cutoff = 0.98;
lz = 120;vol = 120*53*53;
tollen = 0.75;
approx_err = 1; % To compute approximate error in exponentiated Delta F
bstar = 5.0; %Ratio of (bulk-brush_ht)/deb_len


%% Color Scheme

green = [0 0.5 0.0]; gold = [0.9 0.75 0];
orange = [0.91 0.41 0.17]; brown=[0.6 0.2 0];
pclr = {'m',brown,green,'k','b', gold};
lsty = {'-','--',':'};
msty = {'d','s','o','x'};


%% Compute average brush height

fout = fopen('./../../all_txtfiles/htcutoff_all.dat','w');
fprintf(fout,'N\t ht (cutoff = %g)\n',cutoff);
fprintf('N\t ht (cutoff = %g)\n',cutoff);
htvals = zeros(length(nvalsarr),4);

for i = 1:length(nvalsarr)
    
    nval = nvalsarr(i);
    
    for typeval = 1:4
        
        if typeval == 1
            dirstr = 'bl_bl';
        elseif typeval == 2
            dirstr = 'bl_al';
        elseif typeval == 3
            dirstr = 'al_bl';
        else
            dirstr = 'al_al';
        end
        
        fid = fopen(sprintf('./../../results_dens/new_results/results_%d_%s/grpdens.lammpstrj',nval,dirstr));
        data = textscan(fid,'%f%f%f','Headerlines',1);
        fld = cell2mat(data);
        zdata   = fld(:,1);
        pegraft = fld(:,2);
        maxdenval = max(pegraft); cutoffval = (1-cutoff)*maxdenval;
        
        %spline fit
        zspline = 0:0.01:max(zdata);
        denspline = spline(zdata,pegraft,zspline);
        pval = 0;
        for j = 1:length(zspline)-1
            if(denspline(j+1) <= cutoffval && denspline(j) >= cutoffval)
                pval = j;
                break;
            end
        end
        if pval == 0
            disp('Could not find the right height')
        end
        
        ht_cut = 0.5*(zspline(pval)+zspline(pval+1));
        htvals(i,typeval) = ht_cut;
        fprintf(fout,'%g\t%s\t%g\n',nvalsarr(i), dirstr, ht_cut);
        fprintf('%g\t%s\t%g\n',nvalsarr(i), dirstr, ht_cut);
        fclose(fid);
        
    end
    
end
fclose(fout);


h1 = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('$N_{pa}/N_{pc}$','FontSize',20,'Interpreter','Latex')
ylabel('$h_g$ ($\sigma$)','FontSize',20,'Interpreter','Latex')

for i = 1:4
    plot(nvalsarr/ngraft,htvals(:,i),'color',pclr{i},'LineWidth',2, ...
        'LineStyle',lsty{3},'Marker',msty{i},'MarkerSize',8,'MarkerFaceColor',pclr{i})
end

legendinfo{1} = 'Block-Block';
legendinfo{2} = 'Block-Alter';
legendinfo{3} = 'Alter-Block';
legendinfo{4} = 'Alter-Alter';

legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
legend boxoff
saveas(h1,'./../../all_figures/brushht','png');


%% Plot Free Energy

fout = fopen('./../../all_txtfiles/deltaF_all.dat','w');
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
        expfree_energy = zeros(10,3);
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
        
        errbulk = errbulk/sqrt(indfree);
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
saveas(h1,'./../../all_figures/delta_F','png');


%% Exponentially (Boltzmann) averaged Free Energy Plots

free_energy = zeros(10,2);
diff_energy = zeros(length(nvalsarr),4);
err_energy  = zeros(length(nvalsarr),4);
avg_energy  = zeros(length(nvalsarr),4);

fdebye = fopen('./../../all_txtfiles/debye_all.dat','w');
fprintf(fdebye,'%s\t%s\n','Nfree','DebyeLength');

fout = fopen('./../../all_txtfiles/expdeltaF_all.dat','w');
fprintf(fout,'%s\t%s\t%s\t%s\t%s\n','Nfree','Arch','CutoffLength','deltaF','Error');

fht = fopen('./../../all_txtfiles/htcutoff_all.dat','r');
if fht <= 0
    fprintf('%s\n', 'ERROR: No height cutoff file found');
    return;
else
    fgetl(fht);
end

for nvals = 1:length(nvalsarr)
    
    %compute debye length = sqrt(V/(4*pi*(ns+ncounter)))
    
    deb_len = sqrt(vol/(4*pi*(2*nvalsarr(nvals)*nchargedmons+2*ngraft*nchargedmons + nsalt*2)));
    fprintf(fdebye,'%d\t%g\n',nvalsarr(nvals),deb_len);
    
    for i = 1:4
        
        tline = fgetl(fht);
        strarr = strsplit(tline);
        
        ht_cut = str2double(strarr{3}); %brush height
        
        if i == 1
            dirstr = 'bl_bl';
        elseif i == 2
            dirstr = 'bl_al';
        elseif i == 3
            dirstr = 'al_bl';
        else
            dirstr = 'al_al';
        end
        
        fprintf('Boltzmann averaged WHAM analysis for %d,\t %s\n',nvalsarr(nvals),dirstr);
        fylename = sprintf('./../../whamout_all/wham_%d_%s/whamout.txt',nvalsarr(nvals),dirstr);
        fid = fopen(fylename,'r');
        free_energy = zeros(1,3);
        
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
                if strcmp(strarr{3},'-nan') || strcmp(strarr{3},'nan')
                    fprintf('Could not compute errorbar at point %d in %d\t%s\n',k,nvalsarr(nvals),dirstr);
                    fprintf('%s\n','Resetting error to 0');
                    free_energy(k,3) = 0.0;
                else
                    free_energy(k,3) = str2double(strarr{3});
                end
                k = k + 1;
            end
            
        end
        
        [minfree,indmin] = min(free_energy(:,2));
        [maxfree,indmax] = max(free_energy(:,2));
        
        %Check if the 100% of maximum value is less than brushht+bstar*debyelength. If
        %so use the bulk as 95% to 100% or else use bstar*debyelength-100%
        
        brush_lenval = ht_cut + bstar*deb_len; %cutoff length for brush region
        
        lenfree_energy = length(free_energy(:,2));
        flag_totlen = -1;
        if free_energy(lenfree_energy,1) < brush_lenval + bstar*deb_len %bulk should be atleast equal to the width of brush_ht-lencutoff
            flag_totlen = 1;
        end
        
        %Compute bulk free energy and bulk error
        if flag_totlen == 1
            bulkmin_ind = floor(0.95*lenfree_energy);         
        else
            %find nearest value to bstar*debye_len
            for findval = 1:lenfree_energy-1
                if free_energy(findval,1) >= brush_lenval 
                    bulkmin_ind = findval;
                    break;
                end
            end            
        end
        
        bulkfree_energy = mean(free_energy(bulkmin_ind:lenfree_energy,2));
        errbulk = 0.0; indfree = bulkmin_ind; bulkcntr = 0;
        while indfree <= lenfree_energy
            errbulk = free_energy(indfree,3)^2 + errbulk;
            bulkcntr = bulkcntr+1;
            indfree = indfree + 1;
        end
        errbulk = sqrt(errbulk/bulkcntr);
        
        %define deltaF such that deltaF->0 as r->infinity
        %Define new free energy as the difference between standard and
        %bulk
        free_energy(:,2) = free_energy(:,2)-bulkfree_energy;
        expfree_energy = zeros(lenfree_energy,3);
        for j = 1:lenfree_energy
            expfree_energy(j,1) = free_energy(j,1);
            expfree_energy(j,2) = exp(-(free_energy(j,2)));
            expfree_energy(j,3) = sqrt(free_energy(j,3)^2+errbulk^2);
        end
        
        
        %Find the index closest to cutoff length of brush region
        
        flag_ind = -1;
        for lfind = 1:lenfree_energy-1
            if (expfree_energy(lfind,1) >= brush_lenval-tollen) && (expfree_energy(lfind+1,1) < brush_lenval+tollen)
                lcutindex = lfind;
                flag_ind = 1;
                fprintf('lcut at %d\t%g\n',lcutindex,expfree_energy(lcutindex,1));
                break;
            end
        end
        
        if flag_ind == -1
            fprintf('ERROR: Did not find the L0 index for %d\t %s\n ',nvalsarr(nvals),dirstr);
            fprintf('Tolerance length = %g\t, Debye_length = %g\t, brush ht = %g\t, cutoff length = %g\n',...
                tollen,deb_len,ht_cut,brush_lenval);
            break;
        end
        
        %Integrate numerator and denominator. Partition function has to be
        %integrated to the same cutoff length for numerator and denominator
        
        num_val = 0.0;
        den_val = 0.0;
        
        for int_index = 1:lfind
            dx = expfree_energy(int_index+1,1) - expfree_energy(int_index,1);
            dy_num = free_energy(int_index,2)*expfree_energy(int_index,2)+ ...
                free_energy(int_index+1,2)*expfree_energy(int_index+1,2);
            dy_den = expfree_energy(int_index,2)+expfree_energy(int_index+1,2);
            
            num_val = num_val + 0.5*dx*dy_num;
            den_val = den_val + 0.5*dx*dy_den;
        end
              
        avg_energy(nvals,i)  = num_val/den_val;
        diff_energy(nvals,i) = avg_energy(nvals,i);
        
        
        %Computing error.
        %<DeltaF>=sum(a_i b_i)/sum(q_i), where a_i =
        %fz-avgbulk,qi=bi=exp(-ai)
        
        if approx_err == 1
            net_del_of_dden = 0;net_del_of_dnum = 0;
            int_index = 1;
            
            
            dden = expfree_energy(int_index,2);
            del_of_dden = dden*expfree_energy(int_index,3);
            net_del_of_dden = net_del_of_dden + del_of_dden^2;
            
            dnum = free_energy(int_index,2)*expfree_energy(int_index,2);
            net_del_of_dnum  = net_del_of_dnum + dnum^2*(expfree_energy(int_index,3)^2*(1+1/free_energy(int_index,2)^2));
            
            for int_index = 2:lfind-1
                
                dden = 2.0*expfree_energy(int_index,2);
                del_of_dden = dden*expfree_energy(int_index,3);
                net_del_of_dden = net_del_of_dden + del_of_dden^2;
                
                dnum = 2.0*free_energy(int_index,2)*expfree_energy(int_index,2);
                net_del_of_dnum  = net_del_of_dnum + dnum^2*(expfree_energy(int_index,3)^2*(1+1/(free_energy(int_index,2))^2));
                
            end
            
            dnum = free_energy(int_index,2)*expfree_energy(int_index,2);
            net_del_of_dnum  = net_del_of_dnum + dnum^2*(expfree_energy(int_index,3)^2*(1+1/free_energy(int_index,2)^2));
            
            dden = expfree_energy(int_index,2);
            del_of_dden = dden*expfree_energy(int_index,3);
            net_del_of_dden = net_del_of_dden + del_of_dden^2;
            
            err_energy(nvals,i)  = 1/lenfree_energy*(sum(free_energy(:,3).^2));% avg_energy(nvals,i)*sqrt(net_del_of_dden + net_del_of_dnum);
            fprintf(fout,'%d\t%s\t%g\t%g\t%g\n',nvalsarr(nvals),dirstr,brush_lenval,...
                diff_energy(nvals,i),err_energy(nvals,i));
            
        else
            
            err_energy(nvals,i) = sqrt(free_energy(minfree,3)^2 + errbulk);
            fprintf(fout,'%d\t%s\t%g\t%g\t%g\n',nvalsarr(nvals),dirstr,brush_lenval,...
                diff_energy(nvals,i),err_energy(nvals,i));
            
        end
        
    end
    
    fclose(fid);
    
end
fclose(fdebye);
fclose(fout);


h1 = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('$N_{pa}/N_{pc}$','FontSize',20,'Interpreter','Latex')
ylabel('$\langle \Delta F \rangle$ ($k_B T$)','FontSize',20,'Interpreter','Latex')

for i = 1:4
    errorbar(nvalsarr/ngraft,avg_energy(:,i),err_energy(:,i),'color',pclr{i},'LineWidth',2, ...
        'LineStyle',lsty{3},'Marker',msty{i},'MarkerSize',8,'MarkerFaceColor',pclr{i})
end

legendinfo{1} = 'Block-Block';
legendinfo{2} = 'Block-Alter';
legendinfo{3} = 'Alter-Block';
legendinfo{4} = 'Alter-Alter';

legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
legend boxoff
saveas(h1,'./../../all_figures/delta_expF','png');

%% Free_Energy Plots as a fn(z/L)

if delFvsfofz
    
    nvalsarr = [32;48;64;72;80];
    
    for nvals = 1:length(nvalsarr)
        
        h1 = figure;
        hold on
        box on
        set(gca,'FontSize',16)
        xlabel('$z$','FontSize',20,'Interpreter','Latex')
        ylabel('$\Delta F(z)$ ($k_B T$)','FontSize',20,'Interpreter','Latex')
        
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
            
            fylename = sprintf('./whamout_all/new_data/wham_%d_%s/whamout.txt',nvalsarr(nvals),dirstr);
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
            
            errorbar(free_energy(:,1),free_energy(:,2),free_energy(:,3),'color',pclr{i},'LineWidth',2, ...
                'LineStyle',lsty{2},'Marker',msty{i},'MarkerSize',8,'MarkerFaceColor',pclr{i})
            
        end
        
        fclose(fid);
        legendinfo{1} = 'Block-Block';
        legendinfo{2} = 'Block-Alter';
        legendinfo{3} = 'Alter-Block';
        legendinfo{4} = 'Alter-Alter';
        xlim([3 65])
        
        legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','NorthWest')
        legend boxoff
        saveas(h1,sprintf('./../../all_figures/delta_Fz_%d',nvalsarr(nvals)),'png');
        
    end
    
end


%% Histogram Plots

if histflag == 1
    
    nvalsarr = [64;72;80]; 
    
    for nvals = 1:length(nvalsarr)
        
        for i = 1:4
            
            h1 = figure;
            hold on
            box on
            set(gca,'FontSize',16)
            xlabel('$z$','FontSize',20,'Interpreter','Latex')
            ylabel('Probability','FontSize',20,'Interpreter','Latex')
            
            
            if i == 1
                dirstr = 'bl_bl';
            elseif i == 2
                dirstr = 'bl_al';
            elseif i == 3
                dirstr = 'al_bl';
            else
                dirstr = 'al_al';
            end
            
            dirname = sprintf('./whamout_all/new_data/traj_colvars/n_%d_%s',...
                nvalsarr(nvals),dirstr);
            fylepattern = fullfile(dirname,'out.colvars.traj_*'); %Change to reqd pattern.
            dumplist    = dir(fylepattern);
            nfyles = length(dumplist);
            
            for fylecnt = 1:nfyles
                
                baseFileName = dumplist(fylecnt).name;
                fullFileName = fullfile(dirname, baseFileName);
                fprintf('Analyzing %s \n', fullFileName);
                
                fid = fopen(fullFileName,'r');
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
            saveas(h1,sprintf('./../../all_figures/hist_%d_%s',nvalsarr(nvals),dirstr),'png');
            
        end
        
    end

end


%%Plot debye length
h1 = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('$N_{pa}/N_{pc}$','FontSize',20,'Interpreter','Latex')
ylabel('$\lambda_D$ ($\sigma$)','FontSize',20,'Interpreter','Latex')
deb_data = importdata('./../../all_txtfiles/debye_all.dat');
plot(deb_data.data(:,1),deb_data.data(:,2),'--ks','LineWidth',2, ...
    'MarkerSize',8,'MarkerFaceColor','k')

