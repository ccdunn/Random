P_draw = .26;
P_home = .51;
P_away = .23;

g_home = [P_away P_draw 0 P_home];
g_away = [P_home P_draw 0 P_away];

N = 20;

N_home = 17;
N_away = 17;

d_home = g_home;
for ii = 1:N_home
    d_home = conv(d_home,g_home);
end

d_away = g_away;
for ii = 1:N_home
    d_away = conv(d_away,g_away);
end

d_game = conv(d_home,d_away);
plot(0:length(d_game)-1,d_game);

mls_15 = [60 53 51 51 50 49 44 37 37 30 ...
            60 53 53 51 51 51 47 42 41 37];

hold on;
[n,b] = hist(mls_15,[1:length(d_game)]);
stem(b,n./sum(n));
hold off;