function [X1,Y1,Z1,maxk1]=InterpolateGappedTracings(X,Y,Z,maxk,mode,Max_num_dots)

J=length(maxk);
X1=X;
if mode == 0
    Y1=0;
else
    Y1=Y;
end
Z1=Z;
maxk1=maxk;

%calculate the length along each PF
dXt=(diff(X));
dZt=(diff(Z));
if mode == 1
    dYt=(diff(Y));
    dL=sqrt(dXt.^2+dYt.^2+dZt.^2);
else
    dL=sqrt(dXt.^2+dZt.^2);
end
dLfixed=2.5;
flag=0;

%Interpolate PFs with gapped tracing

for j=1:J
    XX=X(1,j);
    if mode == 1
        YY=Y(1,j);
    end
    ZZ=Z(1,j);
    for i=1:maxk(j)-1
        x=[];
        y=[];
        if dL(i,j)>2*dLfixed %spacing between points is more than 2*dLfixed nm
            flag=1;
            points_to_add=round(dL(i,j)/dLfixed);
            x=interp1(X(1:maxk(j),j),i:1/points_to_add:i+1,'pchip');
            if mode==1
                y=interp1(Y(1:maxk(j),j),i:1/points_to_add:i+1,'pchip');
                YY=[YY;y(2:end-1)'];
            end
            z=interp1(Z(1:maxk(j),j),i:1/points_to_add:i+1,'pchip');
            XX=[XX;x(2:end-1)'];
            ZZ=[ZZ;z(2:end-1)'];
        end
        XX=[XX;X(i+1,j)'];
        if mode==1
            YY=[YY;Y(i+1,j)'];
        end
        ZZ=[ZZ;Z(i+1,j)'];
    end
    if flag==1
        maxk1(j)=length(XX);
        X1(1:maxk1(j),j)=XX;
        if mode == 1
            Y1(1:maxk1(j),j)=YY;
        end
        Z1(1:maxk1(j),j)=ZZ;
    end
end
if length(X1(:,1))>Max_num_dots
    X1(Max_num_dots:end,:)=[];
    Z1(Max_num_dots:end,:)=[];
    if mode==1
        Y1(Max_num_dots:end,:)=[];
    end
end
for j=1:J 
    if maxk1(j)>Max_num_dots
        disp(['WARNING: PF#',num2str(j),' was trimmed because longer than ',num2str(Max_num_dots),'points due to interpolatio']);
        maxk1(j)=min(maxk1(j),Max_num_dots);
    end
end


