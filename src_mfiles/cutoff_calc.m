%% To compute brush interface based on some cutoff for density value.

clear;
clc;
close all;
format long;

%% Plot Ion Density Profiles

cutoff = 0.9;
nmons = [32;64;80;100;150]; lz = 120;
typeval = 1; %1 - bl_bl, 2-bl_al, 3-al_bl, 4-al_al

if typeval == 1
    dirstr = 'bl_bl';
elseif typeval == 2
    dirstr = 'bl_al';
elseif typeval == 3
    dirstr = 'al_bl';
else
    dirstr = 'al_al';
end

fout = fopen(sprintf('htcutoff_%s.dat',dirstr),'w');
fprintf(fout,'N\t ht (cutoff = %g)\n',cutoff);
for i = 1:length(nmons)
    
    fid = fopen(sprintf('./densprof/%s/grp_%d.txt',dirstr,nmons(i)));
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
    fprintf(fout,'%g\t%g\n',nmons(i), ht_cut);
    fclose(fid);
    
end
fclose(fout);

% Computing for different architecture
n = 100;lz = 120;
typeval = 1; %1 - bl_bl, 2-bl_al, 3-al_bl, 4-al_al

fout = fopen(sprintf('htcutoff_%d.dat',n),'w');
fprintf(fout,'Arch\t ht (cutoff = %g)\n',cutoff);

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
    
    fid = fopen(sprintf('./densprof/n%d/grp_%d.txt',n,i));
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
    fprintf(fout,'%s\t%g\n',dirstr, ht_cut);
    fclose(fid);
end
fclose(fout);
