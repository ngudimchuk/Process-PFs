%%
% Inputs
folder='.\'; % Results will be saved in this folder

J_dial = "Enter number of protofilaments : ";
J = input(J_dial);

% r_dial = "Enter spacing between tracing points, nm : ";
% r = input(r_dial);
r=2; %minimal spacing between tracing points

c_dial = "Enter curvature, deg/dimer : ";
c = input(c_dial);

sd_dial = "Enter SD of curvature, deg/dimer : ";
sd = input(sd_dial);

curv_grad_dial = "Enter gradient of curvature, deg/dimer^2 : ";
curv_grad = input(curv_grad_dial);

% Array allocation
ang=zeros(100,J);
% c0=-c*pi/180/8*r; % starting curvature, rad/nm
% ang(1,:)=c0;
X=zeros(100,J);
Y=zeros(100,J);
Z=zeros(100,J);
maxk=zeros(1,J);
error_flag=0; 

d_flag = 1; % 1 to delete PFs out of range; 0 to shorten to the boundary

dist_dialog="Select type of PFs length distribution:  fixed / uniform / exponential / normal / gamma / : ";
dist_len=input(dist_dialog, "s"); % Distribution of PF length,
% options include: 'uniform', 'normal', 'exponential', 'gamma', 'fixed'

disp(['Distribution selected: ', dist_len]);

switch dist_len
    case 'uniform'
        lp = "Enter lower PF length limit, nm : ";
        a = input(lp);
        if a<3*r
            disp(['Warning: lower PF length limit set to 3 tracing points =',num2str(3*r),' nm']);
            a=3*r;
        end
        up = "Enter upper PF length limit, nm : ";
        b = input(up);
        if b>98*r
            disp(['Warning: upper PF length limit set to 98 tracing points =',num2str(98*r),' nm']);
            b=98*r;
        end
        len_pf=r*floor(unifrnd(a,b,J,1)/r);
    case 'normal'
        m = "Enter PF length mean value, nm : ";
        a = input(m);
        s = "Enter PF length sigma value, nm : ";
        b = input(s);
        len_pf=floor(normrnd(a,b,J,1)/r)*r;
        outhigh=sum(len_pf>98*r);
        if outhigh>0
            if d_flag >0
                disp(['Warning: N=',num2str(outhigh),' PFs deleted because of exceeding upper length limit (98 tracing points) =',num2str(98*r),' nm']);
                len_pf(len_pf>98*r)=NaN;
            else
                disp(['Warning: N=',num2str(outhigh),' PF lengths set to the upper limit (98 tracing points) =',num2str(98*r),' nm']);
                len_pf(len_pf>98*r)=98*r;
            end
        end
        outlow=sum(len_pf<3*r);
        if outlow>0
            if d_flag>0
                disp(['Warning: N=',num2str(outlow), ' PF deleted, being beyond lower length limit (3 tracing points) =',num2str(3*r),' nm']);
                len_pf(len_pf<3*r)=NaN;
            else
                disp(['Warning: N=',num2str(outlow), ' PF lengths set to the lower limit (3 tracing points) =',num2str(3*r),' nm']);
                len_pf(len_pf<3*r)=3*r;
            end
        end
    case 'exponential'
        m = "Enter mean PF length value, nm : ";
        a = input(m);
        len_pf=floor(exprnd(a,J,1)/r)*r;
        outhigh=sum(len_pf>98*r);
        if outhigh>0
            if d_flag >0
                disp(['Warning: N=',num2str(outhigh),' PFs deleted because of exceeding upper length limit (98 tracing points) =',num2str(98*r),' nm']);
                len_pf(len_pf>98*r)=NaN;
            else
                disp(['Warning: N=',num2str(outhigh),' PF lengths set to the upper limit (98 tracing points) =',num2str(98*r),' nm']);
                len_pf(len_pf>98*r)=98*r;
            end
        end
        outlow=sum(len_pf<3*r);
        if outlow>0
            if d_flag>0
                disp(['Warning: N=',num2str(outlow), ' PF deleted, being beyond lower length limit (3 tracing points) =',num2str(3*r),' nm']);
                len_pf(len_pf<3*r)=NaN;
            else
                disp(['Warning: N=',num2str(outlow), ' PF lengths set to the lower limit (3 tracing points) =',num2str(3*r),' nm']);
                len_pf(len_pf<3*r)=3*r;
            end
        end
    case 'gamma'
        sh = "Enter the shape parameter value for PF length distribution (e.g. 3): ";
        a = input(sh);
        sc = "Enter the scale parameter value for PF length distribution (e.g. 20), nm: ";
        b = input(sc);
        len_pf=floor(gamrnd(a,b,J,1)/r)*r;
        outhigh=sum(len_pf>98*r);
        if outhigh>0
            if d_flag >0
                disp(['Warning: N=',num2str(outhigh),' PFs deleted because of exceeding upper length limit (98 tracing points) =',num2str(98*r),' nm']);
                len_pf(len_pf>98*r)=NaN;
            else
                disp(['Warning: N=',num2str(outhigh),' PF lengths set to the upper limit (98 tracing points) =',num2str(98*r),' nm']);
                len_pf(len_pf>98*r)=98*r;
            end
        end
        outlow=sum(len_pf<3*r);
        if outlow>0
            if d_flag>0
                disp(['Warning: N=',num2str(outlow), ' PF deleted, being beyond lower length limit (3 tracing points) =',num2str(3*r),' nm']);
                len_pf(len_pf<3*r)=NaN;
            else
                disp(['Warning: N=',num2str(outlow), ' PF lengths set to the lower limit (3 tracing points) =',num2str(3*r),' nm']);
                len_pf(len_pf<3*r)=3*r;
            end
        end
    case 'fixed'
        pd = "Enter PF length, nm : ";
        p=input(pd);        
        if p>=98*r
            disp(['Warning: Value exceeds the limit. PF length set to 98 tracing points =',num2str(98*r),' nm']);
            p=98*r;
        end
        if p<3*r
            disp(['Warning: Value exceeds the limit. PF length set to 3 tracing points =',num2str(3*r),' nm']);
            p=3*r;
        end 
        len_pf=repmat(floor(p/r)*r,J,1);
    otherwise
        error_flag=1;
        disp('ERROR: INVALID DISTRIBUTION TYPE SPECIFIED');        
end

if error_flag<1
%%
    for j=1:J
        %add an additional tracing point to mimick MT wall vector in IMOD
        %HOWFLARED routine:
        X_on_MT_wall= 0;
        Y_on_MT_wall= 0;
        Z_on_MT_wall= -5;
        X(1,j)=X_on_MT_wall;
        Y(1,j)=Y_on_MT_wall;
        Z(1,j)=Z_on_MT_wall;
        X(2,j)=0;
        Y(2,j)=0;
        Z(2,j)=0;
        i=2;
        tot_len=0;
        while (round(tot_len/r)*r<=len_pf(j))&&(i<100)
            ang(i,j)=ang(i-1,j)+(c+sd*randn(1))*(pi/180/8)*r+curv_grad*(pi/180/8/8)*(i*r)*r;
            X(i+1,j)=X(i,j)+r*sin(ang(i,j));
            Z(i+1,j)=Z(i,j)+r*cos(ang(i,j));
            tot_len=tot_len+sqrt((X(i+1,j)-X(i,j))^2+(Z(i+1,j)-Z(i,j))^2);
            i=i+1;
        end
        maxk(j)=i;
    end
    
    %bring the first PF tracing point to zero: {X,Y,Z}=(0,0,0)
    for j=1:J
        X(1:maxk(j),j)=X(1:maxk(j),j)-X(1,j);
        Y(1:maxk(j),j)=Y(1:maxk(j),j)-Y(1,j);
        Z(1:maxk(j),j)=Z(1:maxk(j),j)-Z(1,j);
    end
    
    %delete PFs whose lengths were beyond 3 to 98 tracing points limits
    X(:,isnan(len_pf))=[];
    Y(:,isnan(len_pf))=[];
    Z(:,isnan(len_pf))=[];
    maxk(isnan(len_pf))=[];
    J=J-sum(isnan(len_pf));

%     for jj=1:J
%         delpoint=round(rand(1,round(rand(1)*95+1))*98+2);
%         if ~isempty(delpoint)
%             for ii=delpoint
%                 X(ii:end-1,jj)=X(ii+1:end,jj);
%                 Z(ii:end-1,jj)=Z(ii+1:end,jj);
%             end
%         end
%     end
    
    figure
    for j=1:J
        plot(X(1:maxk(j),j),Z(1:maxk(j),j),'k');
        hold on
%         axis square
        axis equal
    end
    
    mkdir([folder,'simulated_coordinates']);
    writetable(table(X),[folder,'simulated_coordinates/x.txt'],'Delimiter','\t','WriteVariableNames',0);
    writetable(table(Z),[folder,'simulated_coordinates/z.txt'],'Delimiter','\t','WriteVariableNames',0);
    writetable(table(Y),[folder,'simulated_coordinates/y.txt'],'Delimiter','\t','WriteVariableNames',0);
end


