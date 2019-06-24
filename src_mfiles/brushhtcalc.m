%% To compute brush height from first moment

clear;
clc;
close all;
format long;

%% Plot Ion Density Profiles

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

fout = fopen(sprintf('brushht_%s.dat',dirstr),'w');
fprintf(fout,'%s\t%s\n','N','h');
for i = 1:length(nmons)
    
    fid = fopen(sprintf('./densprof/%s/grp_%d.txt',dirstr,nmons(i)));
    data = textscan(fid,'%f%f%f','Headerlines',1);
    fld = cell2mat(data);
    zdata   = fld(:,1);
    pegraft = fld(:,2);
    pefree  = fld(:,3);
    
    z_rhog = zdata.*pegraft;
    br_num = trapz(zdata,z_rhog);
    br_den = trapz(zdata,pegraft);
    br_ht  = br_num/br_den;
    
    max(pegraft.*pefree)
    
    fprintf(fout,'%g\t%g\n',nmons(i), br_ht);
    fclose(fid);
    
end
fclose(fout);

% Computing for different architecture
n = 100;lz = 120;
typeval = 1; %1 - bl_bl, 2-bl_al, 3-al_bl, 4-al_al

fout = fopen(sprintf('brushht_%d.dat',n),'w');
fprintf(fout,'%s\t%s\n','Arch','h');
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
    pefree  = fld(:,3);
    
    z_rhog = zdata.*pegraft;
    br_num = trapz(zdata,z_rhog);
    br_den = trapz(zdata,pegraft);
    br_ht  = br_num/br_den;
    
    fprintf(fout,'%s\t%g\n',dirstr, br_ht);
    fclose(fid);
end
fclose(fout);
