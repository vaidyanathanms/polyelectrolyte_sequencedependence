% To compute number of free polymers from density profile

close all;
clear;
clc;
format long;

%% Plot Polymer Density Profiles - ratio of negative to positive charged polymers

% Graft
nchains = 100;cutoff = 0.9;lz = 120;nmons = 30;area=53^2;
pclr = {'r','b','k','m','g'};
nclr = {'--r','--b','--k','--m','--g'};
cclr = {'r*','b*','k*','m*','g*'};

fout = fopen(sprintf('integ_htcutoff_anioncationcharge_rat%d.dat',nchains),'w');
fprintf(fout,'%s\n','Ratio of charged anions to charged cations within cutoff defined by 10% of charged cations (Mikes Method)');
fprintf(fout,'Arch\t ht (cutoff = %g)\t PosIntegral\t NegIntegral\t Rat\n',cutoff);

h4 = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('$r/L_z$','FontSize',20,'Interpreter','Latex')
ylabel('$\rho^{g}(r)$','FontSize',20,'Interpreter','Latex')

for i = 1:4

    if i == 1
        dirstr = 'blbl';
    elseif i == 2
        dirstr = 'blal';
    elseif i == 3
        dirstr = 'albl';
    else
        dirstr = 'alal';
    end
    
    fid = fopen(sprintf('./newresults/dens_%s.txt',dirstr));
    data = textscan(fid,'%f%f%f','Headerlines',1);
    fld = cell2mat(data);
    zdata  = fld(:,1);
    pepos  = fld(:,2);
    peneg  = fld(:,3);
    plot(zdata/lz,pepos, pclr{i}, 'LineWidth', 2, 'MarkerSize', 8)
    
    maxdenval = max(pepos); cutoffval = (1-cutoff)*maxdenval;
    
    %spline fit
    zspline = 0:0.01:max(zdata);
    denspline = spline(zdata,pepos,zspline);
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
    
    %nchk = (53^2)*(7120/1920)*trapz(zdata,pegraft)/1920
    brushspline_pos = spline(zdata,pepos,zspline);
    freespline_neg = spline(zdata,peneg,zspline);
    integpos = 0.0;
    integneg = 0.0;
    for j = 2:pval
        integpos = integpos+0.5*(zspline(j)-zspline(j-1))*(brushspline_pos(j)+brushspline_pos(j+1));
        integneg = integneg+0.5*(zspline(j)-zspline(j-1))*(freespline_neg(j)+freespline_neg(j+1));
    end
         
    integpos = area*integpos*64*15;
    integneg = area*integneg*100*15;
    rat = integneg/integpos;
    ht_cut = 0.5*(zspline(pval)+zspline(pval+1));
    fprintf(fout,'%s\t%g\t%g\t%g\t%g\n',dirstr, ht_cut, integpos, integneg, rat);
    fclose(fid);
    
    
end

legendinfo{1} = 'Block-Block';
legendinfo{2} = 'Block-Alter';
legendinfo{3} = 'Alter-Block';
legendinfo{4} = 'Alter-Alter';

legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
legend boxoff
saveas(h4,'mikemethod','png');

fclose(fout);

%% Plot Polymer Density Profiles - Just integral of polyanions

% Graft
nchains = 100;cutoff = 0.9;lz = 120;
pclr = {'r','b','k','m','g'};

fout = fopen(sprintf('integ_htcutoff_fracanion_%d.dat',nchains),'w');
fprintf(fout,'%s\n','Ratio of total anions to total cations within cutoff defined by 10% of Total cations');
fprintf(fout,'Arch\t ht (cutoff = %g)\t FracPolyAnions\n ',cutoff);

h4 = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('$r/L_z$','FontSize',20,'Interpreter','Latex')
ylabel('$\rho^{g}(r)$','FontSize',20,'Interpreter','Latex')

for i = 1:4

    if i == 1
        dirstr = 'blbl';
    elseif i == 2
        dirstr = 'blal';
    elseif i == 3
        dirstr = 'albl';
    else
        dirstr = 'alal';
    end
    
    fid = fopen(sprintf('./newresults/grp_%s.txt',dirstr));
    data = textscan(fid,'%f%f%f','Headerlines',1);
    fld = cell2mat(data);
    zdata  = fld(:,1);
    pepos  = fld(:,2);
    peneg  = fld(:,3);
    plot(zdata/lz,pepos, pclr{i}, 'LineWidth', 2, 'MarkerSize', 8)
    
    maxdenval = max(pepos); cutoffval = (1-cutoff)*maxdenval;
    
    %spline fit
    zspline = 0:0.01:max(zdata);
    denspline = spline(zdata,pepos,zspline);
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
    
    %nchk = (53^2)*(7120/1920)*trapz(zdata,pegraft)/1920
    brushspline_pos = spline(zdata,pepos,zspline);
    freespline_neg = spline(zdata,peneg,zspline);
    integpos = 0.0;
    integneg = 0.0;
    for j = 2:pval
        integneg = integneg+0.5*(zspline(j)-zspline(j-1))*(freespline_neg(j)+freespline_neg(j+1));
    end
         
    integneg = area*integneg;
    ht_cut = 0.5*(zspline(pval)+zspline(pval+1));
    fprintf(fout,'%s\t%g\t%g\n',dirstr, ht_cut, integneg);
    fclose(fid);
    
    
end

legendinfo{1} = 'Block-Block';
legendinfo{2} = 'Block-Alter';
legendinfo{3} = 'Alter-Block';
legendinfo{4} = 'Alter-Alter';

legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
legend boxoff
saveas(h4,'densgraft','png');
fclose(fout);
%% Plot Polymer Density Profiles 

% Graft
nchains = 100;cutoff = 0.9;lz = 120;
pclr = {'r','b','k','m','g'};

fout = fopen(sprintf('integ_htcutoff_fracchargeanion_%d.dat',nchains),'w');
fprintf(fout,'%s\n','Ratio of Charged anions to total cations within cutoff defined by 10% of Total Cations');
fprintf(fout,'Arch\t ht (cutoff = %g)\t FracPolyAnions\n ',cutoff);

h4 = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('$r/L_z$','FontSize',20,'Interpreter','Latex')
ylabel('$\rho^{g}(r)$','FontSize',20,'Interpreter','Latex')

for i = 1:4

    if i == 1
        dirstr = 'blbl';
    elseif i == 2
        dirstr = 'blal';
    elseif i == 3
        dirstr = 'albl';
    else
        dirstr = 'alal';
    end
    
    fid = fopen(sprintf('./newresults/grp_%s.txt',dirstr));
    data = textscan(fid,'%f%f%f','Headerlines',1);
    fld = cell2mat(data);
    zdata  = fld(:,1);
    pepos  = fld(:,2);
    peneg  = fld(:,3);
    plot(zdata/lz,pepos, pclr{i}, 'LineWidth', 2, 'MarkerSize', 8)
    
    fid2 = fopen(sprintf('./newresults/dens_%s.txt',dirstr));
    data = textscan(fid2,'%f%f%f','Headerlines',1);
    fld = cell2mat(data);
    zdata  = fld(:,1);
    chargepos  = fld(:,2);
    chargeneg  = fld(:,3);
    
    maxdenval = max(pepos); cutoffval = (1-cutoff)*maxdenval;
    
    %spline fit
    zspline = 0:0.01:max(zdata);
    denspline = spline(zdata,pepos,zspline);
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
    
    %nchk = (53^2)*(7120/1920)*trapz(zdata,pegraft)/1920
    brushspline_pos = spline(zdata,chargepos,zspline);
    freespline_neg = spline(zdata,chargeneg,zspline);
    integpos = 0.0;
    integneg = 0.0;
    for j = 2:pval
        integneg = integneg+0.5*(zspline(j)-zspline(j-1))*(freespline_neg(j)+freespline_neg(j+1));
    end
         
    integneg = area*integneg;
    ht_cut = 0.5*(zspline(pval)+zspline(pval+1));
    fprintf(fout,'%s\t%g\t%g\n',dirstr, ht_cut, integneg);
    fclose(fid);fclose(fid2);
    
    
end

legendinfo{1} = 'Block-Block';
legendinfo{2} = 'Block-Alter';
legendinfo{3} = 'Alter-Block';
legendinfo{4} = 'Alter-Alter';

legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
legend boxoff
saveas(h4,'densgraft_2','png');
fclose(fout);

%% Plot Polymer Density Profiles 

% Graft
nchains = 100;cutoff = 0.9;lz = 120;
pclr = {'r','b','k','m','g'};

fout = fopen(sprintf('integ_htcutoff_rat_fracchargeanionfracchargecation_%d.dat',nchains),'w');
fprintf(fout,'%s\n','Ratio of Charged anions to Charged cations within cutoff defined by 10% of Total Cations');
fprintf(fout,'Arch\t ht (cutoff = %g)\t PosIntegral\t NegIntegral\t Rat\n',cutoff);

h4 = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('$r/L_z$','FontSize',20,'Interpreter','Latex')
ylabel('$\rho^{g}(r)$','FontSize',20,'Interpreter','Latex')

for i = 1:4

    if i == 1
        dirstr = 'blbl';
    elseif i == 2
        dirstr = 'blal';
    elseif i == 3
        dirstr = 'albl';
    else
        dirstr = 'alal';
    end
    
    fid = fopen(sprintf('./newresults/grp_%s.txt',dirstr));
    data = textscan(fid,'%f%f%f','Headerlines',1);
    fld = cell2mat(data);
    zdata  = fld(:,1);
    pepos  = fld(:,2);
    peneg  = fld(:,3);
    plot(zdata/lz,pepos, pclr{i}, 'LineWidth', 2, 'MarkerSize', 8)
    
    fid2 = fopen(sprintf('./newresults/dens_%s.txt',dirstr));
    data = textscan(fid2,'%f%f%f','Headerlines',1);
    fld = cell2mat(data);
    zdata  = fld(:,1);
    chargepos  = fld(:,2);
    chargeneg  = fld(:,3);
    
    maxdenval = max(pepos); cutoffval = (1-cutoff)*maxdenval;
    
    %spline fit
    zspline = 0:0.01:max(zdata);
    denspline = spline(zdata,pepos,zspline);
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
    
    %nchk = (53^2)*(7120/1920)*trapz(zdata,pegraft)/1920
    brushspline_pos = spline(zdata,chargepos,zspline);
    freespline_neg = spline(zdata,chargeneg,zspline);
    integpos = 0.0;
    integneg = 0.0;
    for j = 2:pval
        integpos = integpos+0.5*(zspline(j)-zspline(j-1))*(brushspline_pos(j)+brushspline_pos(j+1));
        integneg = integneg+0.5*(zspline(j)-zspline(j-1))*(freespline_neg(j)+freespline_neg(j+1));
    end
         
    integpos = area*integpos;
    integneg = area*integneg;
    rat = integneg/integpos;
    ht_cut = 0.5*(zspline(pval)+zspline(pval+1));
    fprintf(fout,'%s\t%g\t%g\t%g\t%g\n',dirstr, ht_cut, integpos, integneg, rat);
    fclose(fid);fclose(fid2);    
    
end

legendinfo{1} = 'Block-Block';
legendinfo{2} = 'Block-Alter';
legendinfo{3} = 'Alter-Block';
legendinfo{4} = 'Alter-Alter';

legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
legend boxoff
saveas(h4,'densgraft_2','png');
fclose(fout);

%% Plot n(r)



nchains = 100;cutoff = 0.9;lz = 120;
pclr = {'r','b','k','m','g'};

rhofree = nchains*30/(lz*area);


h4 = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('$r$','FontSize',20,'Interpreter','Latex')
ylabel('$n_f(r)$','FontSize',20,'Interpreter','Latex')

for i = 1:4

    if i == 1
        dirstr = 'blbl';
    elseif i == 2
        dirstr = 'blal';
    elseif i == 3
        dirstr = 'albl';
    else
        dirstr = 'alal';
    end
    
    fid = fopen(sprintf('./newresults/PErdf_%s.txt',dirstr));
    data = textscan(fid,'%f%f','Headerlines',1);
    fld = cell2mat(data);
    rdata    = fld(:,1);
    rdfdata  = fld(:,2);
    nofrdata = zeros(length(rdata)-1,1);
    rnrdata = zeros(length(rdata)-1,1);
    rsqgrvals = (rdata.^2).*(rdfdata); 
    
    for j = 1:length(rdata)-1
        nofrdata(j,1) = 0.5*(rdata(j+1)-rdata(j))*(rsqgrvals(j)+rsqgrvals(j+1));
        rnrdata(j,1) = 0.5*(rdata(j+1)+rdata(j));
    end
    
    plot(rnrdata,4*pi*rhofree*nofrdata, pclr{i}, 'LineWidth', 2, 'MarkerSize', 8)
    
    
end

legendinfo{1} = 'Block-Block';
legendinfo{2} = 'Block-Alter';
legendinfo{3} = 'Alter-Block';
legendinfo{4} = 'Alter-Alter';

legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','Best')
legend boxoff
xlim([0 2.5] )
saveas(h4,'PEnofr','png');
fclose(fout);
