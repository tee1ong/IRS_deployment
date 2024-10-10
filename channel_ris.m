function g = channel_ris(nr, nt, nris, imp)
H = randn(nris,nt) + 1j*randn(nris,nt); %tx-ris channel NLOS
h = randn(nr,nris) + 1j*randn(nr,nris); %ris-rx channel
d = randn(nr,nt) + 1j*randn(nr,nt); %tx-rx channel


% phi = ((rand(nris, 1).*abs(imp')) * 2 * pi);
phi = abs(imp')*2*pi;
PHI = diag(exp(1j * phi));
g = (h*PHI*H+d); %overall channel
end