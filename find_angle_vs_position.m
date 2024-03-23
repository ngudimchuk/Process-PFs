function [mean_angle,std_angle, Dist_along_PF, term_ang]=find_angle_vs_position(dA2,dL0,dL,maxk2,Lstep, term_bin)

J=length(maxk2);

for j=1:J
    Lfromwall(1,j)=dL0(j);
    Lfromtip(1,j)=dL0(j);
    for k=2:maxk2(j)
        Lfromwall(k,j)=Lfromwall(k-1,j)+dL(k-1,j);
        Lfromtip(k,j)=Lfromtip(k-1,j)+dL(maxk2(j)-k+1,j);
    end
end
dA2(dA2==0)=NaN;
dAtip=NaN(size(dA2));


for j=1:J
    dAtip(1:maxk2(j)-2,j)=flip(dA2(1:maxk2(j)-2,j));
end

L_edges=2:Lstep:max(max(Lfromtip))+1;
L_bins_fromwall=discretize(Lfromwall,L_edges);
L_bins_fromtip=discretize(Lfromtip,L_edges);

L_bins_fromwall_trimmed=L_bins_fromwall(2:end-1,:);
L_bins_fromtip_trimmed=L_bins_fromtip(2:end-1,:);
L_bins_fromtip_trimmed(isnan(dAtip))=NaN;

term_num=L_bins_fromtip_trimmed==term_bin;  %terminal curvature is determined in the interval from 4 nm to 8 nm (if Lstep=4nm, smoothwindow2=3).
term_ang=mean(dAtip(sum(term_num,2)>0,:),1,'omitnan');

for i=1:length(L_edges-1)
    numwall(i)=sum(sum(L_bins_fromwall_trimmed==i));
    meancurvfromwall(i)=mean(dA2(L_bins_fromwall_trimmed==i),'omitnan');
    stdcurvfromwall(i)=std(dA2(L_bins_fromwall_trimmed==i),'omitnan')./sqrt(numwall(i));
    numtip(i)=sum(sum(L_bins_fromtip_trimmed==i));
    meancurvfromtip(i)=mean(dAtip(L_bins_fromtip_trimmed==i),'omitnan');
    stdcurvfromtip(i)=std(dAtip(L_bins_fromtip_trimmed==i),'omitnan')./sqrt(numtip(i));
end

Dist_along_PF1=L_edges(1:end-1)+Lstep/2;

MinN=2; %minimal required number of PF points to obtain mean value
mean_angle=meancurvfromtip(numtip>MinN);
std_angle=stdcurvfromtip(numtip>MinN);
Dist_along_PF=Dist_along_PF1(numtip>MinN);

