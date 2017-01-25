Ntrials = 1000000;
Nguesses = 10;
A = 1;
B = 1000;
v = randi([A,B],[1 Ntrials]);

r = zeros(1,Ntrials);
for ii = 1:Ntrials
    a = A;
    b = B;
    for jj = 1:Nguesses
        n = Nguesses - jj + 1;
        g = trench(n,a,b);
        
        if g==v(ii)
            r(ii) = v(ii);
        elseif g<v(ii)
            a = g+1;
        else
            b = g-1;
        end
    end
end

mean(r)