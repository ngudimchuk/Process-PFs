function [Xb,Yb,maxkb] = generate_bootstrap_PF_set(X,Y,maxk)

J=length(maxk);
for j=1:J
    rnd=randi(J);
    Xb(:,j)=X(:,rnd);
    Yb(:,j)=Y(:,rnd);
    maxkb(j)=maxk(rnd);
end

