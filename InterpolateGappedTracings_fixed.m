function [X1,Y1,Z1,maxk1]=InterpolateGappedTracings_fixed(X,Y,Z,maxk)

J=length(maxk);
X1=X;
Y1=Y;
Z1=Z;
maxk1=maxk;

%calculate the length along each PF
dXt=(diff(X));
dYt=(diff(Y));
dZt=(diff(Z));
dL=sqrt(dXt.^2+dYt.^2+dZt.^2);
dLfixed=2.5;

%Interpolate PFs with gapped tracing

for j=1:J
    XX=X(1,j);  
    YY=Y(1,j);
    ZZ=Z(1,j);
    for i=1:maxk(j)-1
        if dL(i,j)>2*dLfixed %spacing between points is more than 5nm
            points_to_add=round(dL(i,j)/dLfixed)-1;      
            xstep=dXt(i,j)/(points_to_add+1);
            ystep=dYt(i,j)/(points_to_add+1);
            zstep=dZt(i,j)/(points_to_add+1);
            for k=1:points_to_add
                x(k)=X(i,j)+xstep*k;
                y(k)=Y(i,j)+ystep*k;
                z(k)=Z(i,j)+zstep*k;
            end
            XX=[XX;x'];
            YY=[YY;y'];
            ZZ=[ZZ;z'];
            x=[]; y=[]; z=[];
        end
        XX=[XX;X(i+1,j)'];
        YY=[YY;Y(i+1,j)'];
        ZZ=[ZZ;Z(i+1,j)'];
    end
    maxk1(j)=length(XX);
    X1(1:maxk1(j),j)=XX;
    Y1(1:maxk1(j),j)=YY; 
    Z1(1:maxk1(j),j)=ZZ; 
end



