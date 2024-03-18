function [Xb,Zb,maxkb] = generate_bootstrap_PF_set(X,Y,Z,maxk)

J=length(maxk);
for j=1:J
    rnd=randi(J);
    Xb(:,j)=X(:,rnd);
%     Yb(:,j)=Y(:,rnd);
    Zb(:,j)=Z(:,rnd);
    maxkb(j)=maxk(rnd);
end

