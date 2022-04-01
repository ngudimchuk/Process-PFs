function [dA,dL]=find_angles_sub_test2(Xt,Yt,maxk)
%% Calculating PF lengths

J=length(maxk);
Nz_Yt=Yt~=0;
Ytcum=cumsum(Nz_Yt,1,'reverse');
IndexToNonZero=(Ytcum'>1)';
IndexToNonZero(end,:)=[];

%calculate the length along each PF
dXt=(diff(Xt));
dYt=(diff(Yt));
dL=sqrt(dXt.*dXt+dYt.*dYt);

% remove fragments of zero length (if any)
szL=size(dL);
for j=1:szL(2)
    for k=1:max(maxk)-1
        if (dL(k,j)==0)&&(k<maxk(j)-1)
           [k j] 
           LL=dL(:,j);
           LL=[LL; 0];
           LL(k)=[];
           dL(:,j)=LL;
           XX=dXt(:,j);
           XX=[XX; 0];
           XX(k)=[];
           dXt(:,j)=XX;
           YY=dYt(:,j);
           YY=[YY; 0];
           YY(k)=[];
           dYt(:,j)=YY;
           maxk(j)=maxk(j)-1;
        end
    end
end

%% calculate mean curvature of MT wall proximal PF part
alfa=zeros(max(maxk)-1,J);
for j=1:J
    for k=1:maxk(j)-1
        if (dXt(k,j)>0)&&(dYt(k,j)>0)
            alfa(k,j)=pi/2-acos(abs(dXt(k,j))/dL(k,j));
        elseif (dXt(k,j)>0)&&(dYt(k,j)<=0)
            alfa(k,j)=pi/2+acos(abs(dXt(k,j))/dL(k,j));
        elseif (dXt(k,j)<=0)&&(dYt(k,j)<=0)
            alfa(k,j)=3*pi/2-acos(abs(dXt(k,j))/dL(k,j));
        elseif (dXt(k,j)<=0)&&(dYt(k,j)>0)
            alfa(k,j)=-pi/2+acos(abs(dXt(k,j))/dL(k,j));
        end
    end
end


mL=(dL(1:max(maxk)-2,:)+dL(2:max(maxk)-1,:))/2;
dA=diff(alfa); %rad per interval









