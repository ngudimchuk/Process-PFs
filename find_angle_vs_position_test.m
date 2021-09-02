function [mean_angle, std_angle,Dist_along_PF]=find_angle_vs_position_test(dA2,dL,maxk2)

J=length(maxk2);

for j=1:J
    Lfromwall(1,j)=0;
    Lfromtip(1,j)=0;
    for k=2:maxk2(j)
        Lfromwall(k,j)=Lfromwall(k-1,j)+dL(k-1,j);
        Lfromtip(k,j)=Lfromtip(k-1,j)+dL(maxk2(j)-k+1,j);
    end
end
dAtip=0*dA2;
for j=1:J
    dAtip(1:maxk2(j)-2,j)=flip(dA2(1:maxk2(j)-2,j));
end
Lstep=4;
L_edges=Lstep/2:Lstep:max(max(Lfromtip));
L_bins_fromwall=discretize(Lfromwall,L_edges);
L_bins_fromtip=discretize(Lfromtip,L_edges);

L_bins_fromwall_trimmed=L_bins_fromwall(2:end-1,:);
L_bins_fromtip_trimmed=L_bins_fromtip(2:end-1,:);
for i=1:length(L_edges-1)
    numwall(i)=sum(sum(L_bins_fromwall_trimmed==i));
    meancurvfromwall(i)=nanmean(dA2(L_bins_fromwall_trimmed==i));
    stdcurvfromwall(i)=std(dA2(L_bins_fromwall_trimmed==i),'omitnan')./sqrt(numwall(i));
    numtip(i)=sum(sum(L_bins_fromtip_trimmed==i));
    meancurvfromtip(i)=nanmean(dAtip(L_bins_fromtip_trimmed==i));
    stdcurvfromtip(i)=std(dAtip(L_bins_fromtip_trimmed==i),'omitnan')./sqrt(numtip(i));
end

Dist_along_PF1=L_edges(1:end-1)+Lstep/2;

MinN=2; %minimal required number of PF points to obtain mean value
mean_angle=meancurvfromtip(numtip>MinN);
std_angle=stdcurvfromtip(numtip>MinN);
Dist_along_PF=Dist_along_PF1(numtip>MinN);
% mean_angle=meancurvfromwall(numtip>MinN);
% std_angle=stdcurvfromwall(numtip>MinN);
% Dist_along_PF=Dist_along_PF1(numtip>MinN);