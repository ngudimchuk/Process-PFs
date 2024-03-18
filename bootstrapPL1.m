function [PL,xX,zZ,logsem,f]=bootstrapPL1(Xt,Yt,Zt,maxk,plotflag,bootstrapflag,Lmed,pix,mode)
%
% Xt=Xt1;
% Yt=Yt1;
% maxk=maxk1;
% bootstrapflag=0;
% plotflag=1;
%%

global nameofdataset

if bootstrapflag==1
    [Xt3,Zt3,maxk3] = generate_bootstrap_PF_set(Xt,Yt,Zt,maxk);
else
    Xt3=Xt;
    Yt3=Yt;
    Zt3=Zt;
    maxk3=maxk;
end

[dA3,dL3]=find_angles_sub_test2(Xt3,Zt3,maxk3,mode); %rad per interval


J=length(maxk3);
Lfromtip=zeros(max(maxk3),J);
for j=1:J
    Lfromtip(1,j)=0;
    for k=2:maxk3(j)
        Lfromtip(k,j)=Lfromtip(k-1,j)+dL3(maxk3(j)-k+1,j);
    end
end

dA3(dA3==0)=NaN;
dAtip=NaN(size(dA3));

for j=1:J
    dAtip(1:maxk3(j)-2,j)=flip(dA3(1:maxk3(j)-2,j)); 
end

np=10;
Le=sum(sum(~isnan(dAtip),2)>np);
while isempty(Le)
    np=np/2;
    Le=sum(sum(~isnan(dAtip),2)>np);
end
    
Lstep=Lmed*pix;
L_edges=[Lstep/10:Lstep:Le];%32;
L_bins_fromtip=discretize(Lfromtip,L_edges);
L_bins_fromtip_trimmed=L_bins_fromtip(2:end-1,:);
for i=1:length(L_edges-1)
    numtip(i)=sum(sum(L_bins_fromtip_trimmed==i,'omitnan'),'omitnan');
    meancurvfromtip(i)=mean(dAtip(L_bins_fromtip_trimmed==i),'omitnan');
end

for j=1:J
    alfafromtip(:,j)=cumsum(dAtip(:,j),'omitnan');
end
meanalfatip=cumsum(meancurvfromtip,'omitnan');

for i=1:length(L_edges-1)
    if numtip(i)>10
        delta_ang=(alfafromtip(L_bins_fromtip_trimmed==i)-meanalfatip(i));
        meancosa(i)=mean(cos(delta_ang),'omitnan');
        semcosa(i)=std(cos(delta_ang),'omitnan')./sqrt(numtip(i));
    else
        meancosa(i)=-2;
        semcosa(i)=-2;
    end
end
Dist_along_PF=L_edges(1:end-1)+Lstep/2;

logcos=log(meancosa(meancosa>0));
logsem=semcosa(meancosa>0)./meancosa(meancosa>0);

Weights=1./logsem.^2;
Weights(Weights==Inf)=10^-5;

xX=Dist_along_PF(meancosa>0);
zZ=logcos;
[f,gof] = fit(xX',zZ','-(x-x0)/pl','Weights',Weights,'StartPoint', [100, 0]);
% [f,gof] = fit(xX',yY','-(x-x0)/pl','StartPoint', [1, 100]);

PL=f.pl;

if plotflag==1
    subplot(1,1,1)
    errorbar(xX,zZ,logsem,'b.-');
    hold on
    plot(xX,f(xX),'r-');

    hT1=title([nameofdataset,' - aligned by tip']);
    hx1=xlabel('distance along PF,nm');
    hy1=ylabel('log(<cos(teta)>)');
    hT1.FontSize=18;
    hx1.FontSize=16;
    hy1.FontSize=16;
end
