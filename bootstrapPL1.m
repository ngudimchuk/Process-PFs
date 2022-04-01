function [PL,xX,yY,logsem,f]=bootstrapPL1(Xt, Yt, maxk,plotflag,bootstrapflag,Lstep)
%
% Xt=Xt1;
% Yt=Yt1;
% maxk=maxk1;
% bootstrapflag=0;
% plotflag=1;
%%

global nameofdataset

if bootstrapflag==1
    [Xt3,Yt3,maxk3] = generate_bootstrap_PF_set(Xt,Yt,maxk);
else
    Xt3=Xt;
    Yt3=Yt;
    maxk3=maxk;
end

[dA3,dL3]=find_angles_sub_test2(Xt3,Yt3,maxk3); %rad per interval


J=length(maxk3);
Nz_Yt=Yt3~=0;
Ytcum=cumsum(Nz_Yt,1,'reverse');
IndexToNonZero=(Ytcum'>1)';
IndexToNonZero(end,:)=[];
Lfromtip=zeros(max(maxk3),J);
for j=1:J
    Lfromtip(1,j)=0;
    for k=2:maxk3(j)
        Lfromtip(k,j)=Lfromtip(k-1,j)+dL3(maxk3(j)-k+1,j);
    end
end
dAtip=0*dA3;
for j=1:J
    dAtip(1:maxk3(j)-2,j)=flip(dA3(1:maxk3(j)-2,j)); 
end

L_edges=[Lstep/2:Lstep:32];%max(max(Lfromtip))];
L_bins_fromtip=discretize(Lfromtip,L_edges);
L_bins_fromtip_trimmed=L_bins_fromtip(2:end-1,:);
for i=1:length(L_edges-1)
    numtip(i)=sum(sum(L_bins_fromtip_trimmed==i));
    meancurvfromtip(i)=mean(dAtip(L_bins_fromtip_trimmed==i));
end

for j=1:J
    alfafromtip(:,j)=cumsum(dAtip(:,j));
end
meanalfatip=cumsum(meancurvfromtip);

for i=1:length(L_edges-1)
    if numtip(i)>19
        delta_ang=(alfafromtip(L_bins_fromtip_trimmed==i)-meanalfatip(i));
        meancosa(i)=mean(cos(delta_ang));
        semcosa(i)=std(cos(delta_ang))./sqrt(numtip(i));
    else
        meancosa(i)=-2;
        semcosa(i)=-2;
    end
end
Dist_along_PF=L_edges(1:end-1)+Lstep/2;

logcos=log(meancosa(meancosa>0));
logsem=semcosa(meancosa>0)./meancosa(meancosa>0);

Weights=1./logsem.^2;

xX=Dist_along_PF(meancosa>0);
yY=logcos;
[f,gof] = fit(xX',yY','-(x-x0)/pl','Weights',Weights,'StartPoint', [100, 0]);
% [f,gof] = fit(xX',yY','-(x-x0)/pl','StartPoint', [1, 100]);

PL=f.pl;

if plotflag==1
    subplot(1,1,1)
    errorbar(xX,yY,logsem,'b.-');
    hold on
    plot(xX,f(xX),'r-');
    ax=axis;
    axis([0 ax(2) ax(3) 0]);
    hT1=title([nameofdataset,' - aligned by tip']);
    hx1=xlabel('distance along PF,nm');
    hy1=ylabel('log(<cos(teta)>)');
    hT1.FontSize=18;
    hx1.FontSize=16;
    hy1.FontSize=16;
end