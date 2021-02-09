function [KR,VR,rep] = var2riss(VARA,V,y,r)

ny   = length(y);
nr   = length(r);
p    = size(VARA,3);
pny  = p*ny;
pny1 = pny-ny;

A = [reshape(VARA(y,y,:),ny,pny); eye(pny1) zeros(pny1,ny)];
C = reshape(VARA(r,y,:),nr,pny);
Q = [V(y,y) zeros(ny,pny1); zeros(pny1,pny)];
S = [V(y,r); zeros(pny1,nr)];
R = V(r,r);
[KR,VR,rep] = ss2iss(A,C,Q,R,S);
