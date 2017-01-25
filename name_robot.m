function [name, count] = name_robot(count)

if ~exist('count','var')
    count = 1;
end

P = 2^31 - 1;
Q = 676009;

ind = mod(P*count,Q)

if ind>676000
    [name, count] = name_robot(count + 1);
else
   ind1 = round(ind/(26*1000))
   indtmp = ind - ind1*26000;
   ind2 = round((ind - ind1*26000)/(1000))
   indtmp = indtmp - ind1*26000;
   ind3 = round((ind - ind1*26000 - ind2*1000)/(100))
   indtmp = indtmp - ind1*26000;
   ind4 = round((ind - ind1*26000 - ind2*1000 - ind3*100)/(10))
   indtmp = indtmp - ind1*26000;
   ind5 = round((ind - ind1*26000 - ind2*1000 - ind3*100 - ind4*10)/(1))
   indtmp = indtmp - ind1*26000;
   
   a1 = char(ind1 + 65);
   a2 = char(ind2 + 65);
   
   name = [a1 a2 num2str(ind3) num2str(ind4) num2str(ind5)];
end
   
count = count+1;

end