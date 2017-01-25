Ns = 100;
Nb = 32;

ns = 1:Ns;
ss = 0:Nb-1;

P = zeros(Ns,Nb-1);

[nn,ss] = meshgrid(ns,ss);

P = 2.^(nn.*(ss-Nb)).*(1 - 2.^(1 - nn)).^ss.*factorial(Nb)./(factorial(Nb-ss).*factorial(ss));

E = ns.*sum(P,1);
fprintf('%.10f\n',round(cumsum(E)*10^10)/10^10);
% fprintf('%.10f\n',nansum(P(:)));