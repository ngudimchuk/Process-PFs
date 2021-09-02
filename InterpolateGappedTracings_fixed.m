function [X1,Y1,Z1,maxk1]=InterpolateGappedTracings_fixed(X,Y,Z,maxk,mode)

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

%Interpolate PFs with gapped tracing

for j=1:J
    XX=X(1,j);
    if mode == 1
        YY=Y(1,j);
    end
    ZZ=Z(1,j);
    for i=1:maxk(j)-1
        if dL(i,j)>2*dLfixed %spacing between points is more than 5nm
            points_to_add=round(dL(i,j)/dLfixed)-1;
            xstep=dXt(i,j)/(points_to_add+1);
            if mode == 1
                ystep=dYt(i,j)/(points_to_add+1);
            end
            zstep=dZt(i,j)/(points_to_add+1);
            for k=1:points_to_add
                x(k)=X(i,j)+xstep*k;
                if mode == 1
                    y(k)=Y(i,j)+ystep*k;
                end
                z(k)=Z(i,j)+zstep*k;
            end
            XX=[XX;x'];
            if mode == 1
                YY=[YY;y'];
            end
            ZZ=[ZZ;z'];
            x=[]; y=[]; z=[];
        end
        XX=[XX;X(i+1,j)'];
        if mode == 1
            YY=[YY;Y(i+1,j)'];
        end
        ZZ=[ZZ;Z(i+1,j)'];
    end
    maxk1(j)=length(XX);
    X1(1:maxk1(j),j)=XX;
    if mode == 1
        Y1(1:maxk1(j),j)=YY;
    end
    Z1(1:maxk1(j),j)=ZZ;
end



