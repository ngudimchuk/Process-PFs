function [Xt, Yt, Zt, maxk]=filter_PFs(Xt, Yt, Zt, maxk, smoothwind,stype, addwall, rtips)
J=length(maxk);

% add arificial MT-wall part
if addwall>0
    wallsize=addwall;
    Xwall=0.0001*ones(wallsize,length(maxk));
    Ywall=0.0001*ones(wallsize,length(maxk));
    Zwall=[];
    for i=1:wallsize
        Zwall=[-2.5*i*ones(1,length(maxk)); Zwall];
    end

    Xt=[Xwall;  Xt];
    Yt=[Ywall;  Yt];
    Zt=[Zwall;  Zt];
    maxk=maxk+wallsize;
end

%smooth
smoothdata=stype;
if smoothdata==1 % use simple moving average filter
    for j=1:J
        XX=Xt(1:maxk(j),j);
        Xt(1:maxk(j),j)=smooth(XX,smoothwind);
        
        YY=Yt(1:maxk(j),j);
        Yt(1:maxk(j),j)=smooth(YY,smoothwind);
        
        ZZ=Zt(1:maxk(j),j);
        Zt(1:maxk(j),j)=smooth(ZZ,smoothwind);
    end
elseif smoothdata==2 % use quadratic loess smooth
    for j=1:J
        XX=Xt(1:maxk(j),j);
        Xt(1:maxk(j),j)=smooth(XX,smoothwind/maxk(j),'loess');
        
        YY=Yt(1:maxk(j),j);
        Yt(1:maxk(j),j)=smooth(YY,smoothwind/maxk(j),'loess');
        
        ZZ=Zt(1:maxk(j),j);
        Zt(1:maxk(j),j)=smooth(ZZ,smoothwind/maxk(j),'loess');
    end
end

%remove wall
if addwall>0
    Xt(1:wallsize,:)=[];
    Yt(1:wallsize,:)=[];
    Zt(1:wallsize,:)=[];
    maxk=maxk-wallsize;
end

%remove tips that were not smoothed
removetips=rtips;
if removetips==1
    tipsize=(smoothwind-1)/2;
    for j=1:J
        num=(tipsize+1:maxk(j)-tipsize);
        
        XX=Xt(num,j);
        Xt(:,j)=0.*Xt(:,j);
        Xt(1:length(num),j)=XX;
        
        YY=Yt(num,j);
        Yt(:,j)=0.*Yt(:,j);
        Yt(1:length(num),j)=YY;
        
        ZZ=Zt(num,j);
        Zt(:,j)=0.*Zt(:,j);
        Zt(1:length(num),j)=ZZ;
        
        maxk(j)=length(num);
    end
end



