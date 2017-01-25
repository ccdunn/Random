function [s,ss] = scorefun(truth,entry)

ss = truth - entry;
s = sum(abs(ss));

end