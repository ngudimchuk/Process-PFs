function [dA,dL,dB]=find_angles_sub(Xt,Yt,Zt,maxk,mode)
%% Calculating PF lengths

J=length(maxk);


%calculate the length along each PF
dXt=(diff(Xt));
if mode == 1
    dYt=(diff(Yt));
end
dZt=(diff(Zt));
if mode == 1
    dL=sqrt(dXt.^2+dYt.^2+dZt.^2);
else
    dL=sqrt(dXt.^2+dZt.^2);
end

% remove fragments of zero length (if any)
szL=size(dL);
for j=1:szL(2)
    for k=1:max(maxk)-1
        if (dL(k,j)==0)&&(k<maxk(j)-1)
            LL=dL(:,j);
            LL=[LL; 0];
            LL(k)=[];
            dL(:,j)=LL;
            
            XX=dXt(:,j);
            XX=[XX; 0];
            XX(k)=[];
            dXt(:,j)=XX;
            if mode == 1
                YY=dYt(:,j);
                YY=[YY; 0];
                YY(k)=[];
                dYt(:,j)=YY;
            end
            ZZ=dZt(:,j);
            ZZ=[ZZ; 0];
            ZZ(k)=[];
            dZt(:,j)=ZZ;
            
            maxk(j)=maxk(j)-1;
        end
    end
end

%% calculate mean curvature of MT wall proximal PF part
alfa=NaN(max(maxk)-2,J);
for j=1:J
    for k=2:maxk(j)-1
        if mode == 1
            alfa(k-1,j)=2*((Yt(k,j)-Yt(k-1,j))*(Zt(k+1,j)-Zt(k-1,j))-(Zt(k,j)-Zt(k-1,j))*(Yt(k+1,j)-Yt(k-1,j))+...
                (Zt(k,j)-Zt(k-1,j))*(Xt(k+1,j)-Xt(k-1,j))-(Xt(k,j)-Xt(k-1,j))*(Zt(k+1,j)-Zt(k-1,j))+...
                (Xt(k,j)-Xt(k-1,j))*(Yt(k+1,j)-Yt(k-1,j))-(Yt(k,j)-Yt(k-1,j))*(Xt(k+1,j)-Xt(k-1,j)))/...
                sqrt(((Xt(k,j)-Xt(k-1,j))^2+(Yt(k,j)-Yt(k-1,j))^2+((Zt(k,j)-Zt(k-1,j))^2)) *...
                ((Xt(k+1,j)-Xt(k,j))^2+(Yt(k+1,j)-Yt(k,j))^2+((Zt(k+1,j)-Zt(k,j))^2)) *...
                ((Xt(k-1,j)-Xt(k+1,j))^2+(Yt(k-1,j)-Yt(k+1,j))^2+((Zt(k-1,j)-Zt(k+1,j))^2)));
        else
            alfa(k-1,j)=-2*((Xt(k,j)-Xt(k-1,j))*(Zt(k+1,j)-Zt(k,j)) - (Zt(k,j)-Zt(k-1,j))*(Xt(k+1,j)-Xt(k,j)))/...
                sqrt(((Xt(k,j)-Xt(k-1,j))^2+((Zt(k,j)-Zt(k-1,j))^2)) * ((Xt(k+1,j)-Xt(k,j))^2+((Zt(k+1,j)-Zt(k,j))^2))...
                * ((Xt(k-1,j)-Xt(k+1,j))^2+((Zt(k-1,j)-Zt(k+1,j))^2)));
        end
    end
end

dA=alfa*8/pi*180; %deg per dimer 

%% calculate helix angle
if mode == 1
    beta=zeros(max(maxk)-1,J);
    for j=1:J
        for k=1:maxk(j)-2
            N1 = cross([dXt(k,j),dYt(k,j),dZt(k,j)],[dXt(k+1,j),dYt(k+1,j),dZt(k+1,j)]);
            N2 = cross([dXt(k+1,j),dYt(k+1,j),dZt(k+1,j)],[dXt(k+2,j),dYt(k+2,j),dZt(k+2,j)]);
            beta(k,j)=acosd(round((sum(N1.*N2))/(norm(N1)*norm(N2)),4));
        end
    end
    mL=(dL(1:max(maxk)-2,:)+dL(2:max(maxk)-1,:))/2;
    dB=beta(1:end-1,:)./mL; %deg per nm
else
    dB=0;
end