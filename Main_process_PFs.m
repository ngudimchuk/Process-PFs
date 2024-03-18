%version 3.a March17,2024
function Main_process_PFs(nameofdataset,folder,savefolder,pix,excl_pf,VisQ,view_PFs,mode,term_curv)

Max_num_dots=100; %maximum allowed number of points per PF tracing
if isempty(savefolder)
    mkdir(strcat(folder,'/Results'));
    savefolder=strcat(folder,'/Results');
end


%% LOAD data
% OPTION1: start from the coordinates provided in x.txt, y.txt, z.txt files
if exist([folder,'/x.txt'])&& exist([folder,'/z.txt'])
    Xt_raw=load([folder,'/x.txt']);
    Zt_raw=load([folder,'/z.txt']);
    if exist([folder,'/y.txt'])
        Yt_raw=load([folder,'/y.txt']);
    end
    disp('DATA LOADED FROM USER-DEFINED .txt files')    
% OPTION2:Converting raw PF tracings data from iMOD in .dat format    
else
    
    fnames = dir(strcat(folder,'/*.dat'));
    
    X=[];
    Z=[];
    if mode == 1
        Y=[];
    end
    global_MT_num=[];
    MT_prev_max(1)=0;

    ii=1;
    while isempty(strfind(fnames(ii).name,'start')) %for each POINT file in folder
    
        % read data from POINT file
        fid=fopen([strcat(folder,'/'),fnames(ii).name]);
        fnames(ii).name;
        tline = fgetl(fid);
        tlines = cell(0,1);
        while ischar(tline)
            tlines{end+1,1} = tline;
            tline = fgetl(fid);
        end
        fclose(fid);
        clear info_lines;
        for i=1:length(tlines)
            info_lines(i)=length(str2num(tlines{i}))>3;
        end
        PFinfo1=tlines(info_lines); % cell array
        PFinfo=str2num(cell2mat(PFinfo1)); % regular array
        Coordinates1=tlines(~info_lines); % cell array
        Coordinates=str2num(cell2mat(Coordinates1)); % regular array
    
        % record data from START and POINT files to matrices
        k=0;
        PFnum(ii)=length(PFinfo(:,1));
        if ii>1
            MT_prev_max(ii)=max(global_MT_num);
        end
        MT_nums=transpose(PFinfo(:,5))+MT_prev_max(ii);
        global_MT_num=[global_MT_num MT_nums];
        x=zeros(Max_num_dots,PFnum(ii));
        z=zeros(Max_num_dots,PFnum(ii));
        if mode == 1
            y=zeros(Max_num_dots,PFnum(ii));
        end
    
        for j=1:PFnum(ii)
            num=k+1:k+PFinfo(j,1);
            k=k+PFinfo(j,1);
            x(1:PFinfo(j,1),j)=Coordinates(num,1);
            z(1:PFinfo(j,1),j)=Coordinates(num,2);
            if mode == 1
                y(1:PFinfo(j,1),j)=Coordinates(num,3);
            end
        end
    
        X=[X x];
        Z=[Z z];
    
        if mode == 1
            Y=[Y y];
        end
    
        if length(fnames)>ii
            ii=ii+1;
        else
            break;
        end
    end

    % create a folder to store x,y,z coordinates
    mkdir([folder,'_processed_coordinates']);
    writetable(table(X),[folder,'_processed_coordinates/x.txt'],'Delimiter','\t','WriteVariableNames',0);
    writetable(table(Z),[folder,'_processed_coordinates/z.txt'],'Delimiter','\t','WriteVariableNames',0);
    if mode == 1
        writetable(table(Y),[folder,'_processed_coordinates/y.txt'],'Delimiter','\t','WriteVariableNames',0);
    end


    Xt_raw=X;
    Zt_raw=Z;
    if mode == 1
        Yt_raw=Y;
    end
    disp('DATA LOADED from HOWFLARED OUTPUT FILES (.dat)') 
end
  
%% PROCESS LOADED PF coordinates

% exclude user-specified 'bad' PF tracings (if expludePFs.txt provided)
if exist([folder,'/excludePFs.txt'])
    jexcl2=load([folder,'/excludePFs.txt']);
    disp(['N = ',num2str(length(jexcl2)),' PFs',' were excluded from analysis by expludePFs.txt file']);
else
    jexcl2=[];
end
Xt_raw(:,jexcl2)=[];
if mode==1
    Yt_raw(:,jexcl2)=[];
end
Zt_raw(:,jexcl2)=[];

%delete empty columns, i.e. remove void PF tracings (if any)
jexl=[];
sz=size(Xt_raw);
if mode == 0
    for j=1:sz(2) 
        if sum(Xt_raw(:,j))+sum(Zt_raw(:,j))==0 
            jexl=[jexl j];
        end
    end
else
    for j=1:sz(2) 
        if sum(Xt_raw(:,j))+sum(Yt_raw(:,j))+sum(Zt_raw(:,j))==0 
            jexl=[jexl j];
        end
    end
end
Xt_raw(:,jexl)=[];
Zt_raw(:,jexl)=[];
if mode == 1
    Yt_raw(:,jexl)=[];
end

%determine the number of non-zero coordinates (maxk0) for each PF
if mode == 1
    Nz=(abs(Xt_raw)+abs(Yt_raw)+abs(Zt_raw))~=0;
else
    Nz=(abs(Xt_raw)+abs(Zt_raw))~=0;
end
maxk0=sum(Nz,1);
J=length(maxk0); %number of protofilaments in the dataset


% orient PF to be in the 1st quadrant
Xt0=zeros(Max_num_dots,J);
Yt0=zeros(Max_num_dots,J);
Zt0=zeros(Max_num_dots,J);
for j=1:J
    hypo=(sqrt(Xt_raw(1:maxk0(j)+1,j).^2+Zt_raw(1:maxk0(j)+1,j).^2));
    dh=diff(hypo);
    if maxk0(j)>2
        if abs(dh(1))>abs(mean(dh)+3*std(dh))
            if mode == 1
                Yt_raw(1:min([maxk0(j)+1,Max_num_dots-1]),j)=Yt_raw(2:min([maxk0(j)+2,Max_num_dots]),j);
                Yt_raw(1:maxk0(j),j)=Yt_raw(1:maxk0(j),j)-Yt_raw(1,j);
            end
            Xt_raw(1:min([maxk0(j)+1,Max_num_dots-1]),j)=Xt_raw(2:min([maxk0(j)+2,Max_num_dots]),j);
            Zt_raw(1:min([maxk0(j)+1,Max_num_dots-1]),j)=Zt_raw(2:min([maxk0(j)+2,Max_num_dots]),j);
            Xt_raw(1:maxk0(j),j)=Xt_raw(1:maxk0(j),j)-Xt_raw(1,j);
            Zt_raw(1:maxk0(j),j)=Zt_raw(1:maxk0(j),j)-Zt_raw(1,j);
        end
        [~,inx]=findpeaks(Xt_raw(1:maxk0(j)+1,j),'MinPeakHeight',15); %[~,inx,~,prx]
        [~,inz]=findpeaks(Zt_raw(1:maxk0(j)+1,j),'MinPeakHeight',15); %[~,inz,~,prz]
        [~,inxn]=findpeaks(-Xt_raw(1:maxk0(j)+1,j),'MinPeakHeight',15); %[~,inxn,~,prxn]
        [~,inzn]=findpeaks(-Zt_raw(1:maxk0(j)+1,j),'MinPeakHeight',15); %[~,inzn,~,przn]
        warning('off')
        if isempty(inx) && isempty(inxn)
            inx=maxk0(j)+1;
        else
            inx=[inx;inxn];
            inx=sort(inx);
            inx=inx(1);
        end
        if isempty(inz) && isempty(inzn)
            inz=maxk0(j)+1;
        else
            inz=[inz;inzn];
            inz=sort(inz);
            inz=inz(1);
        end
    else
        inx=maxk0(j)+1;
        inz=maxk0(j)+1;
    end
    
    XCOM=mean(Xt_raw(2:inx,j));
    ZCOM=mean(Zt_raw(2:inz,j));
    Xt0(:,j)=Xt_raw(:,j)*sign(XCOM);
    Zt0(:,j)=Zt_raw(:,j)*sign(ZCOM);
    if mode == 1
        Yt0(:,j)=Yt_raw(:,j)*sign(XCOM)*sign(ZCOM);
    end
end 


%delete the first tracing point (the point used only for the purpose of MT wall vector
%determination in HOWFLARED IMOD program)
Xt0(1,:)=[];
Zt0(1,:)=[];
if mode == 1
   Yt0(1,:)=[];
else
   Yt0=[];
end

%update the number of points (maxk0) for each PF 
if mode == 1
    Nz=(abs(Xt0)+abs(Yt0)+abs(Zt0))~=0;
else
    Nz=(abs(Xt0)+abs(Zt0))~=0;
end
maxk0=sum(Nz,1);

% [Nz(:,431) Xt0(:,431) Zt0(:,431)]

%fill in tracing gaps by interpolation
[Xt0,Yt0,Zt0,maxk]=InterpolateGappedTracings(Xt0,Yt0,Zt0,maxk0,mode,Max_num_dots);

%bring the first point to origin
if mode == 1
    for j=1:length(maxk)
        Xt0(1:maxk(j),j)=Xt0(1:maxk(j),j)-Xt0(1,j);
        Yt0(1:maxk(j),j)=Yt0(1:maxk(j),j)-Yt0(1,j);
        Zt0(1:maxk(j),j)=Zt0(1:maxk(j),j)-Zt0(1,j);
    end
else
    for j=1:length(maxk)
        Xt0(1:maxk(j),j)=Xt0(1:maxk(j),j)-Xt0(1,j);
        Zt0(1:maxk(j),j)=Zt0(1:maxk(j),j)-Zt0(1,j);
    end
end

% [Nz(:,431) Xt0(:,431) Zt0(:,431)]

% [Nz(:,431) Xt0(:,431) Zt0(:,431)]
% Express PF coordinates in nm using user-supplied pix/nm parameter
Xt=pix*Xt0;
if mode == 1
    Yt=pix*Yt0;
else
    Yt=[];
end
Zt=pix*Zt0;


%% CALCULATE PF LENGTHS
[L_init, Lmed]=PF_lengths_sub(Xt, Yt, Zt, maxk, mode);

% OPTIONAL: exclude PFs, which are shorter than a user-defined length
if excl_pf>0  
    jexcl3=find(L_init<excl_pf*8.0);
    
    Xt(:,jexcl3)=[];
    if mode == 1
        Yt(:,jexcl3)=[];
    end
    Zt(:,jexcl3)=[];
    maxk(jexcl3)=[];
    L=L_init;
    L(:,jexcl3)=[];
    disp(['N = ',num2str(length(jexcl3)),' PFs shorter than ',num2str(excl_pf),' dimers were excluded from analysis !']);
else
    L=L_init;
end

fig_len=figure;
lenbins=(0:(max(L)+10)/20:max(L)+10);
histogram(L, lenbins);
[Nlen,len_edges]=histcounts(L,lenbins);
len_bin_middles=(len_edges(1:end-1)+len_edges(2:end))/2;
hT2=title(['PF lengths: ',nameofdataset]);
xlabel('PF length,nm');
ylabel('Number of PFs');
hT2.FontSize=16;
hx2.FontSize=14;
hy2.FontSize=14;
set(gca,'fontsize',14)

saveas(fig_len,[savefolder,'/histlength.jpeg']);
dlmwrite([savefolder,'/all_PF_lengths.txt'],L','delimiter','\t','newline','pc');
dlmwrite([savefolder,'/hist_PF_length.txt'],[len_bin_middles' Nlen'],'delimiter','\t','newline','pc');
dlmwrite([savefolder,'/MeanL_StdL_numL.txt'], [mean(L) std(L) length(L)],'delimiter','\t','newline','pc');


%% plot all PFs before smoothing

figure
for j=1:length(maxk)
%     if (Zt(maxk(j),j)==0)&&(Xt(maxk(j),j)==0)
%         j
%     end
    if mode == 1
        plot3(Xt(1:maxk(j),j),Yt(1:maxk(j),j),Zt(1:maxk(j),j),'k-');
    else
        plot(Xt(1:maxk(j),j),Zt(1:maxk(j),j),'k-');
    end
    hold on
end
axis square
get(gca,'XLim');
hT2=title(nameofdataset);
if mode == 1
    view(-30,20)
    hx=xlabel('x, nm');
    hy=ylabel('y, nm');
    hz=zlabel('distance along MT axis, nm');
    hx.FontSize=16;
    hy.FontSize=16;
    hz.FontSize=14;
    hT2.FontSize=14;
    axis([-10 max(max(Xt))+20 -10 max(max(Yt))+20 -10 max(max(Zt))+20])
else
    hx=xlabel('x, nm');
    hz=ylabel('distance along MT axis, nm');
    hx.FontSize=16;
    hz.FontSize=16;
    hT2.FontSize=14;
    axis([-50 250 -50 250]);
end

set(gca,'fontsize',14)
get(0,'Factory');
set(0,'defaultfigurecolor',[1 1 1]);

%% OPTIONAL: mark some PFs to identify poorly processed PF tracings (available only in 2D)
% bad tracings can be manually added to excludePFs.txt if needed). 

if (VisQ==1)&&(mode==0)
    disp('1)LEFT MOUSE CLICK - to select a PF. 2)LEFT CLICK - to accept, RIGHT - to deselect. 3)PRESS SCROLL - to finish.');
    axis([-10 100 -10 100]);
    d=100*ones(1,length(maxk));
    w=0;
    while w<1
          [x,z,button]=ginput(1);
           if button==1
               for j=1:length(maxk)
                   d(j)=min((x-Xt(1:maxk(j),j)).^2+(z-Zt(1:maxk(j),j)).^2);
               end
               [~,mind_j]=min(d);
               plot(Xt(1:maxk(mind_j),mind_j), Zt(1:maxk(mind_j),mind_j),'r-');
               [~,~,button2]=ginput(1);
               if button2==3
                  plot(Xt(1:maxk(mind_j),mind_j),Zt(1:maxk(mind_j),mind_j),'k-');
               elseif button2==2
                   break;
               else                   
                  disp(['PF# ',num2str(mind_j),' marked']);
               end
           elseif button==2
               break;
           end 
    end
end
                
%% OPTIONAL: view individual PFs before smoothing
if view_PFs==1
    figure
    j=1;
    while j<=length(maxk)
        if mode == 1
            plot3(Xt(1:maxk(j),j),Yt(1:maxk(j),j),Zt(1:maxk(j),j),'ko-');
        else
            plot(Xt(1:maxk(j),j),Zt(1:maxk(j),j),'ko-');
        end
        hold off
        axis square
        grid minor;
        if mode == 1
            NumTicks = 9;
            scL = get(gca,'XLim');
            set(gca,'XTick',linspace(scL(1),scL(2),NumTicks))
            set(gca,'YTick',linspace(scL(1),scL(2),NumTicks))
        end
        title(['PF # ',num2str(j)]);
        axis equal;
        axis square
        grid minor;
        % to go back, press any key. to go forward, mouse click
        w=waitforbuttonpress;
        j=j+1;
        if (w==1)&&(j>2)
           j=j-2;
        end
    end
end

%% SMOOTH PFs
%Smoothing parameters:
smoothwind=10;
addwall=smoothwind;
cuttips=0;
[Xt1, Yt1, Zt1, maxk1]=filter_PFs(Xt, Yt, Zt, maxk, smoothwind,2,addwall,cuttips,mode);
J=length(maxk1);

%bring first point to origin (smoothing process above can slightly shift it)
if mode == 1
    for j=1:length(maxk1)
        Xt1(1:maxk1(j),j)=Xt1(1:maxk1(j),j)-Xt1(1,j);
        Yt1(1:maxk1(j),j)=Yt1(1:maxk1(j),j)-Yt1(1,j);
        Zt1(1:maxk1(j),j)=Zt1(1:maxk1(j),j)-Zt1(1,j);
    end
else
    for j=1:length(maxk1)
        Xt1(1:maxk1(j),j)=Xt1(1:maxk1(j),j)-Xt1(1,j);
        Zt1(1:maxk1(j),j)=Zt1(1:maxk1(j),j)-Zt1(1,j);
    end
end

% plot all PFs after smoothing
fig_PFs2=figure;
for j=1:length(maxk1)
    if mode == 1
        plot3(Xt1(1:maxk1(j),j),Yt1(1:maxk1(j),j),Zt1(1:maxk1(j),j),'k-');
    else
        plot(Xt1(1:maxk1(j),j),Zt1(1:maxk1(j),j),'k-');
    end
    hold on
end
axis square
hT2=title(nameofdataset);
if mode == 1
    view(-30,20)
    hx=xlabel('x, nm');
    hy=ylabel('y, nm');
    hz=zlabel('distance along MT axis, nm');
    hx.FontSize=16;
    hy.FontSize=16;
    hz.FontSize=14;
    axis([-10 max(max(Xt1))+20 -10 max(max(Yt1))+20 -10 max(max(Zt1))+20]);
else
    hx=xlabel('x, nm');
    hz=ylabel('distance along MT axis, nm');
    hx.FontSize=16;
    hz.FontSize=16;
    axis([-50 250 -50 250]);
end

hT2.FontSize=14;
set(gca,'fontsize',14)
saveas(fig_PFs2,[savefolder,'/plot_PFs.jpeg']);


%% CALCULATE PF angles
% apply smoothing
smoothwind=10;
addwall=smoothwind;
cuttips=0;
[Xt1, Yt1, Zt1, maxk1]=filter_PFs(Xt, Yt, Zt, maxk, smoothwind,2,addwall,cuttips,mode);
J=length(maxk1);

[dA,~,~]=find_angles_sub(Xt1,Yt1,Zt1,maxk1,mode); %deg/dimer

%pool all angles across all PFs
dA_all=[];
for j=1:J
    dA_all=[dA_all (dA(1:maxk1(j)-2,j))'];
end

curvbins=(max([-90,(min(dA_all)-10)]):min([90,max(dA_all)+10]));
[N_dA_all,curv_edges]=histcounts(dA_all,curvbins);

curv_bin_middles=(curv_edges(1:end-1)+curv_edges(2:end))/2;

fig_ang=figure;

histogram(dA_all,curvbins);
hT1=title(['All angles: ', nameofdataset]);
hx1=xlabel('Angle,deg/dimer');
hy1=ylabel('Occurence');
hT1.FontSize=14;
hx1.FontSize=14;
hy1.FontSize=14;
set(gca,'fontsize',14)

dlmwrite([savefolder,'/all_angles.txt'],dA_all','delimiter','\t','newline','pc');
dlmwrite([savefolder,'/Ang-ALL_mean_med_SD_N.txt'],[mean(dA_all(abs(dA_all)<180))...
    median(dA_all(abs(dA_all)<180)) std(dA_all(abs(dA_all)<180))...
    length(dA_all(abs(dA_all)<180))],'delimiter','\t','newline','pc');
dlmwrite([savefolder,'/hist_ang-all.txt'],[curv_bin_middles' N_dA_all'],'delimiter','\t','newline','pc');

saveas(fig_ang,[savefolder,'/angles.jpeg']);

%% mean angle as a function of position along PF
%smooth data with moving average filter with trimmed ends (milder smoothing
%is usually needed)
smoothwind2=3;
addwall=0;
cuttips=1;
[Xt2, Yt2, Zt2, maxk2]=filter_PFs(Xt, Yt, Zt, maxk, smoothwind2,1,addwall,cuttips,mode);
[dA2,dL,~]=find_angles_sub(Xt2,Yt2,Zt2,maxk2,mode); %deg/dimer

% overlay smoothed and cut PFs used for curvature vs distance analysis upon
% smoothed uncut PFs
plot_flag=0;
if plot_flag==1
    figure;
    for j=1:length(maxk1)
        if mode == 1
            plot3(Xt1(1:maxk1(j),j),Yt1(1:maxk1(j),j),Zt1(1:maxk1(j),j),'k-');
            hold on
            plot3(Xt2(1:maxk2(j),j),Yt2(1:maxk2(j),j),Zt2(1:maxk2(j),j),'b-');
        else
            plot(Xt1(1:maxk1(j),j),Zt1(1:maxk1(j),j),'k-');
            hold on
            plot(Xt2(1:maxk2(j),j),Zt2(1:maxk2(j),j),'b-');
        end
        
    end
end

%shift along PF by dL0 is the point at the tip was cut
dL0=zeros(1,length(maxk2));
if cuttips==1
    tipsize=(smoothwind2-1)/2;
    dXt=diff(Xt);
    if mode==1
        dYt=diff(Yt);
    end
    dZt=diff(Xt);
    for i=1:tipsize
        if mode==1
            dL0=dL0+sqrt(dXt(i,:).^2+dYt(i,:).^2+dZt(i,:).^2);
        else
            dL0=dL0+sqrt(dXt(i,:).^2+dZt(i,:).^2);
        end
    end
else
    dL0=zeros(1,length(maxk2));
end

Lstep=4; %nm, binning of the curvature vs distance curve
term_bin=1; % the bin, where terminal curvature is collected. In case of Lstep=4nm, 
% 'term_bin=2' means that the 'terminal curvature' is defined in the (4nm; 8nm] range 

[meancurvfromtip, stdcurvfromtip, Dist_along_PF, term_ang]=find_angle_vs_position(dA2,dL0,dL,maxk2,Lstep, term_bin);
err=stdcurvfromtip(1:length(Dist_along_PF));

Xfit=Dist_along_PF;
Zfit=meancurvfromtip(1:length(Dist_along_PF));
Efit=err;
Xfit(isnan(err))=[];
Zfit(isnan(err))=[];
Efit(isnan(err))=[];

res=Xfit<=80;
Xfit=Xfit(res);
Zfit=Zfit(res);
Efit=Efit(res);

Efit(Efit==0)=0.05;

[f,~] = fit(Xfit',Zfit','sl*x+x0','Weights',1./Efit.^2,'StartPoint',[1 1]);

xnum=[0:1:80];
figure
errorbar(Dist_along_PF,meancurvfromtip(1:length(Dist_along_PF)),stdcurvfromtip(1:length(Dist_along_PF)));
h=errorbar(Xfit,Zfit,Efit);
hold on
plot(xnum,f(xnum),'k-','LineWidth',1);

h.LineWidth=2;
h.Color=[1 0 0];
h.Marker='.';
h.MarkerSize=5;
grid on
axis([0 80 0 40]); 
hT1=title(['Aligned by tip: ', nameofdataset]);
hx1=xlabel('Distance from tip, nm');
hy1=ylabel('Angle,deg/dimer');
hT1.FontSize=14;
hx1.FontSize=14;
hy1.FontSize=14;
set(gca,'fontsize',14)

dlmwrite([savefolder,'/ang_vs_tip.txt'],[Dist_along_PF',...
    (meancurvfromtip(1:length(Dist_along_PF)))',...
    (stdcurvfromtip(1:length(Dist_along_PF)))'],'delimiter','\t','newline','pc');
dlmwrite([savefolder,'/ang_vs_tip_fit.txt'],[xnum' (f(xnum))],'delimiter','\t','newline','pc');

if term_curv>0
    disp(['mean terminal curvature: ', num2str(round(meancurvfromtip(term_bin),1)),' +/- ', num2str(round(stdcurvfromtip(1),1)), ' deg/dimer']);
    fig_termang=figure;
    dlmwrite([savefolder,'/term_angles.txt'],(term_ang(~isnan(term_ang)))','delimiter','\t','newline','pc');
    histogram(term_ang(~isnan(term_ang)),curvbins);
    hT1=title(['terminal angles: ', nameofdataset]);
    hx1=xlabel('terminal angle,deg/dimer');
    hy1=ylabel('Occurence');
    hT1.FontSize=14;
    hx1.FontSize=14;
    hy1.FontSize=14;
    set(gca,'fontsize',14)
    saveas(fig_termang,[savefolder,'/term_angles.jpeg']);
end
