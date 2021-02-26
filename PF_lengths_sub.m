function [L, dLmedian]=PF_lengths_sub(Xt, Yt, Zt, maxk)

%% Calculating PF lengths
J=length(maxk);
Nz=(abs(Xt)+abs(Yt)+abs(Zt))~=0;
Ytcum=cumsum(Nz,1,'reverse');
IndexToNonZero=(Ytcum'>1)';
IndexToNonZero(end,:)=[];

%calculate the length along each PF
dXt=(diff(Xt));
dYt=(diff(Yt));
dZt=(diff(Zt));
dL=sqrt(dXt.^2+dYt.^2+dZt.^2);

dLmedian=median(median(dL(IndexToNonZero)));

%remove fragments of zero length
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
           
           YY=dYt(:,j);
           YY=[YY; 0];
           YY(k)=[];
           dYt(:,j)=YY;
           
           ZZ=dZt(:,j);
           ZZ=[ZZ; 0];
           ZZ(k)=[];
           dZt(:,j)=ZZ;
           
           maxk(j)=maxk(j)-1;
        end
    end
end

L=max((cumsum(dL.*IndexToNonZero,1)));