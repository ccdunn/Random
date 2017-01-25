N = 4;

N/12*sin(pi/N)*(cos(pi/N) + 2)/(1 + cos(pi/N))^2


th = 0:.001:pi/4;
r = 1./(1 + cos(th));

[x,y] = pol2cart(th,r);

plot(x,y,'k');
hold on;
plot(x,-y,'k');
plot(-x,-y,'k');
plot(-x,y,'k');
plot(y,x,'k');
plot(y,-x,'k');
plot(-y,-x,'k');
plot(-y,x,'k');
plot([1 1 -1 -1 1],[1 -1 -1 1 1],'k');
hold off;
aet;