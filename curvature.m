function [L,R,k] = curvature(X)
  N = size(X,1);
  dims = size(X,2);
  if dims == 2
    X = [X,zeros(N,1)];
  end
  L = zeros(N,1);
  R = NaN(N,1);
  k = NaN(N,3);
  D = NaN(N,3);
  for i = 2:N-1
    [R(i),~,k(i,:),D(i,:)] = circ(X(i,:)',X(i-1,:)',X(i+1,:)');
    L(i) = L(i-1)+norm(X(i,:)-X(i-1,:));
  end
  i = N;
  L(i) = L(i-1)+norm(X(i,:)-X(i-1,:));
  if dims == 2
    k = k(:,1:2);
  end
end