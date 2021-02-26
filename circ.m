function [R,M,k,D] = circ(a1,a2,a3)
  D = cross(a2-a1,a3-a1);
  b = norm(a1-a3);
  c = norm(a1-a2);
  if nargout == 1
    a = norm(a2-a3); 
    R = a*b*c/2/norm(D);
    return
  end
  E = cross(D,a2-a1);
  F = cross(D,a3-a1);  
  G = (b^2*E-c^2*F)/norm(D)^2/2;
  M = a1 + G;
  R = norm(G); 
  if R == 0
    k = G;
  else
    k = G'/R^2;   
  end
end