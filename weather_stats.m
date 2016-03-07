temps = [92	92	91	89	86	85	83	82	81	80	79;...
78	79	79	78	77	75	71	69	67	65	64;...
71	72	73	72	71	68	65	62	61	59	58;...
65	66	66	65	64	62	60	58	56	55	54;...
58	59	58	58	58	56	54	52	50	49	48;...
53	53	53	53	52	51	50	49	47	46	45;...
44	43	44	42	40	40	39	39	40	38	38];

times = .5 + [2:12]/24 + datenum('April 23 2016');

fill([times fliplr(times) times(1)],[temps(7,:) fliplr(temps(1,:)) temps(7,1)],[.99 .99 .99],'EdgeColor','none');

hold on;
h(9) = fill([times fliplr(times) times(1)],[temps(6,:) fliplr(temps(2,:)) temps(6,1)],[.9 .9 .9],'EdgeColor','none');
h(8) = fill([times fliplr(times) times(1)],[temps(5,:) fliplr(temps(3,:)) temps(5,1)],[.7 .7 .7],'EdgeColor','none');

linex = [times(1); times(end)];
lineys = 40:10:100;
line(repmat(linex,1,numel(lineys)),repmat(lineys,2,1),'LineWidth',1,'Color',[.85 .85 .85]);
linexs = times;
liney = [32; 100];
line(repmat(linexs,2,1),repmat(liney,1,numel(linexs)),'LineWidth',1,'Color',[.85 .85 .85]);
h(7) = plot(times,temps(1,:),'Color',[1 .9 .9]);
h(6) = plot(times,temps(2,:),'Color',[.95 .7 .7]);
h(5) = plot(times,temps(3,:),'Color',[.9 .4 .4]);
h(4) = plot(times,temps(4,:),'Color','k');
h(3) = plot(times,temps(5,:),'Color',[.4 .4 .9]);
h(2) = plot(times,temps(6,:),'Color',[.6 .7 .95]);
h(1) = plot(times,temps(7,:),'Color',[.9 .9 1]);
hold off;

axis([min(times) max(times) 32 100]);

datetick('x','HH PM');
set(gca,'position',[.1 .1 .8 .8]);
grid off;
GridInFront = 1;

legend(h,{'Minimum','10% Percentile','25% Percentile','Mean','75% Percentile',...
    '90% Percentile','Maximum','Central 50% range','Central 80% range'});
title('Historical Weather Data for April 23 near Ellicott City, MD');
xlabel('Time of Day on April 23');
ylabel([char(176) 'F']);