function err = ladder(x)

l1 = 20;
l2 = 30;
h = 10;

h1 = sqrt(l1.^2 - x.^2);
h2 = sqrt(l2.^2 - x.^2);

xe = h.*( sqrt(l1.^2./h1.^2 - 1) + sqrt(l2.^2./h2.^2 - 1) );

err = abs(x - xe);

end