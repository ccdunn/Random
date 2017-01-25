
clear all;
import_data;

ngames = (size(nycfc2016,1) - 2)/2;

stats.earliestgf = intmax;
stats.latestgf = 0;
stats.earliestga = intmax;
stats.latestga = 0;

stats.scorenext = zeros(8,8,3);
stats.scoreresult = zeros(8,8,3);
stats.scoreends = zeros(8,8);
stats.gas = [];
stats.gfs = [];
stats.gainds = [];
stats.gfinds = [];

for ii = 1:ngames
    
    game(ii).date = datetime(nycfc2016{2*ii+1,1},'InputFormat','MM/dd/yy');
    game(ii).against = nycfc2016{2*ii+2,2};
    
    game(ii).home = strcmp(nycfc2016{2*ii+2,1},'HOME');
    game(ii).away = ~game(ii).home;
    
    game(ii).gf = nycfc2016{2*ii+1,3};
    game(ii).ga = nycfc2016{2*ii+2,3};
    game(ii).gd = game(ii).gf - game(ii).ga;
    
    game(ii).win = game(ii).gf>game(ii).ga;
    game(ii).draw = game(ii).gf==game(ii).ga;
    game(ii).loss = game(ii).gf<game(ii).ga;
    game(ii).points = 3*game(ii).win + 1*game(ii).draw + 0*game(ii).loss;
    
    jj = 7;
    game(ii).gfs = zeros(1,game(ii).gf);
    while ~isempty(nycfc2016{2*ii+1,jj})
        tmp = nycfc2016{2*ii+1,jj};
        if ischar(tmp)
            if strfind(tmp,'+')
                game(ii).gfs(jj-6) = [1 .1]*[sscanf(tmp,'%i+%i')];
            else
                game(ii).gfs(jj-6) = sscanf(tmp,'%i');
            end
        else
            game(ii).gfs(jj-6) = tmp;
        end
        jj = jj + 1;
    end
    jj = 7;
    game(ii).gas = zeros(1,game(ii).ga);
    while ~isempty(nycfc2016{2*ii+2,jj})
        tmp = nycfc2016{2*ii+2,jj};
        if ischar(tmp)
            if strfind(tmp,'+')
                game(ii).gas(jj-6) = [1 .1]*[sscanf(tmp,'%i+%i')];
            else
                game(ii).gas(jj-6) = sscanf(tmp,'%i');
            end
        else
            game(ii).gas(jj-6) = tmp;
        end
        jj = jj + 1;
    end
    
    game(ii).losing = 0;
    game(ii).winningconcede = 0;
    game(ii).drawingconcede = 0;
    game(ii).losingconcede = 0;
    for jj = 1:length(game(ii).gas)
        if game(ii).gas(jj)<stats.earliestga
            stats.earliestga = game(ii).gas(jj);
        end
        if game(ii).gas(jj)>stats.latestga
            stats.latestga = game(ii).gas(jj);
        end
        if jj>sum(game(ii).gfs<game(ii).gas(jj))
            game(ii).losing = 1;
        end
        if (jj-1)>sum(game(ii).gfs<game(ii).gas(jj))
            game(ii).winningconcede = game(ii).winningconcede + 1;
        elseif (jj-1)==sum(game(ii).gfs<game(ii).gas(jj))
            game(ii).drawingconcede = game(ii).drawingconcede + 1;
        elseif (jj-1)<sum(game(ii).gfs<game(ii).gas(jj))
            game(ii).losingconcede = game(ii).losingconcede + 1;
        end
    end
    
    game(ii).winning = 0;
    game(ii).winningscore = 0;
    game(ii).drawingscore = 0;
    game(ii).losingscore = 0;
    for jj = 1:length(game(ii).gfs)
        if game(ii).gfs(jj)<stats.earliestgf
            stats.earliestgf = game(ii).gfs(jj);
        end
        if game(ii).gfs(jj)>stats.latestgf
            stats.latestgf = game(ii).gfs(jj);
        end
        if jj>sum(game(ii).gas<game(ii).gfs(jj))
            game(ii).winning = 1;
        end
        if (jj-1)>sum(game(ii).gas<game(ii).gfs(jj))
            game(ii).winningscore = game(ii).winningscore + 1;
        elseif (jj-1)==sum(game(ii).gas<game(ii).gfs(jj))
            game(ii).drawingscore = game(ii).drawingscore + 1;
        elseif (jj-1)<sum(game(ii).gas<game(ii).gfs(jj))
            game(ii).losingscore = game(ii).losingscore + 1;
        end
    end
    
    gfjj = 1;
    gajj = 1;
    score = [0 0];
    resind = 1*game(ii).win + 2*game(ii).draw + 3*game(ii).loss;
    stats.scoreresult(score(1)+1,score(2)+1,resind) = stats.scoreresult(score(1)+1,score(2)+1,resind) + 1;
    for jj = 1:game(ii).gf + game(ii).ga
        if gajj>game(ii).ga || (gfjj<=game(ii).gf && game(ii).gfs(gfjj)<game(ii).gas(gajj))
            stats.scorenext(score(1)+1,score(2)+1,1) = stats.scorenext(score(1)+1,score(2)+1,1) + 1;
            score = score + [1 0];
            gfjj = gfjj + 1;
        elseif  gfjj>game(ii).gf || (gajj<=game(ii).ga && game(ii).gfs(gfjj)>game(ii).gas(gajj))
            stats.scorenext(score(1)+1,score(2)+1,2) = stats.scorenext(score(1)+1,score(2)+1,2) + 1;
            score = score + [0 1];
            gajj = gajj + 1;
        else
            keyboard;
        end
        stats.scoreresult(score(1)+1,score(2)+1,resind) = stats.scoreresult(score(1)+1,score(2)+1,resind) + 1;
    end
    stats.scorenext(score(1)+1,score(2)+1,3) = stats.scorenext(score(1)+1,score(2)+1,3) + 1;
    
    game(ii).comebackwin = game(ii).losing && game(ii).win;
    game(ii).comebackdraw = game(ii).losing && game(ii).draw;
    game(ii).aheadloss = game(ii).winning && game(ii).loss;
    game(ii).aheaddraw = game(ii).winning && game(ii).draw;
    
    game(ii).wescorefirst = game(ii).gf>0 && (game(ii).ga==0 || (game(ii).gfs(1) < game(ii).gas(1)));
    game(ii).theyscorefirst = game(ii).ga>0 && (game(ii).gf==0 || (game(ii).gas(1) < game(ii).gfs(1)));
    
    if ii==1
        game(ii).dayssincelastplayed = datenum(game(ii).date - datetime('10/25/15','InputFormat','MM/dd/yy'));
    else
        game(ii).dayssincelastplayed = datenum(game(ii).date - game(ii-1).date);
    end
    
    stats.scoreends(game(ii).gf+1,game(ii).ga+1) = stats.scoreends(game(ii).gf+1,game(ii).ga+1) + 1;
    
    game(ii).firstgoal = min(cat(1,game(ii).gfs(:),game(ii).gas(:),100));
    if isempty(game(ii).firstgoal)
        game(ii).firstgoal = nan;
    end
    
    stats.gas = cat(1,stats.gas,game(ii).gas(:));
    stats.gainds = cat(1,stats.gainds,ii*ones(length(game(ii).gas(:)),1));
    
    stats.gfs = cat(1,stats.gfs,game(ii).gfs(:));
    stats.gfinds = cat(1,stats.gfinds,ii*ones(length(game(ii).gfs(:)),1));
    
end

stats.scorepredict = (stats.scoreresult(:,:,:))./repmat(sum(stats.scoreresult,3),[1 1 3]);

stats.gdaxis = -7:7;
for ii = stats.gdaxis
    for jj = 1:3
        stats.gdpredict(ii+8,jj) = sum(diag(stats.scoreresult(:,:,jj),-ii))./sum(diag(sum(stats.scoreresult,3),-ii));
    end
end
imagesc(1:3,stats.gdaxis,stats.gdpredict);
aet;
set(gca,'YDir','normal');


stats.points = sum([game(:).points]);
stats.gf = sum([game(:).gf]);
stats.ga = sum([game(:).ga]);
stats.gd = sum([game(:).gd]);

stats.winpercentofwinning = sum([game(logical([game(:).winning])).win])./sum([game(:).winning]);
stats.drawpercentofwinning = sum([game(logical([game(:).winning])).draw])./sum([game(:).winning]);
stats.losspercentofwinning = sum([game(logical([game(:).winning])).loss])./sum([game(:).winning]);
stats.winpercentoflosing = sum([game(logical([game(:).losing])).win])./sum([game(:).losing]);
stats.drawpercentoflosing = sum([game(logical([game(:).losing])).draw])./sum([game(:).losing]);
stats.losspercentoflosing = sum([game(logical([game(:).losing])).loss])./sum([game(:).losing]);

stats.winpercentofwescorefirst = sum([game(logical([game(:).wescorefirst])).win])./sum([game(:).wescorefirst]);
stats.drawpercentofwescorefirst = sum([game(logical([game(:).wescorefirst])).draw])./sum([game(:).wescorefirst]);
stats.losspercentofwescorefirst = sum([game(logical([game(:).wescorefirst])).loss])./sum([game(:).wescorefirst]);
stats.winpercentoftheyscorefirst = sum([game(logical([game(:).theyscorefirst])).win])./sum([game(:).theyscorefirst]);
stats.drawpercentoftheyscorefirst = sum([game(logical([game(:).theyscorefirst])).draw])./sum([game(:).theyscorefirst]);
stats.losspercentoftheyscorefirst = sum([game(logical([game(:).theyscorefirst])).loss])./sum([game(:).theyscorefirst]);


pearson = @(x,y) mean((x - mean(x)).*(y - mean(y)))./sqrt(var(x).*var(y));

fns = fieldnames(game(1));
isnum = ones(1,length(fns));
for ii = 1:length(fns)
    for jj = 1:ngames
        if ~isnumeric(game(jj).(fns{ii})) && ~islogical(game(jj).(fns{ii})) || length(game(jj).(fns{ii}))>1
            isnum(ii) = 0;
            break;
        end
    end
end
isnum = find(isnum);

corr = nan(length(isnum));
for ii = 1:length(isnum)
    for jj = ii+1:length(isnum)
        corr(ii,jj) = pearson([game(:).(fns{isnum(ii)})],[game(:).(fns{isnum(jj)})]);
    end
end

[scorr,scind] = sort(abs(colvec(corr)),'descend');
scind = scind(~isnan(scorr));
scorr = scorr(~isnan(scorr));

[scindx, scindy] = ind2sub([length(isnum),length(isnum)],scind);

corrs = cat(2,fns(isnum(scindx)),fns(isnum(scindy)),num2cell(corr(scind)));
corrs

% http://sport.maths.org/content/ball-0
% So in the Premiership, indeed most professional soccer, we expect a team to win about 2/3 of the games in which it scores first, and draw about 1/5 of them. That offers the warm comfort that if your team scores first, it should lose only about one time in seven. You can check the match outcomes each week, and over a season, from information in the newspapers. Real data do conform well to these
% proportions.

% figure(1);
% plot([game.date],cumsum([game.points]));
% hold on;
% plot([game.date],cumsum([game.gf]));
% plot([game.date],cumsum([game.ga]));
% hold off;
% print('./output/SeasonByDate.png','-dpng');
% 
% figure(2);
% plot(cumsum([game.points]));
% hold on;
% plot(cumsum([game.gf]));
% plot(cumsum([game.ga]));
% hold off;
% legend({'points','goals for','goals against'});
% print('./output/SeasonByGame.png','-dpng');
% 
% figure(3);
% imagesc([0:7],[0:7],squeeze(stats.scorepredict(:,:,1)));
% aet;
% set(gca,'XAxisLocation','top');
% ylabel('NYCFC Goals');
% xlabel('Opponent Goals');
% title('Win Probability');
% print('./output/WinProbability.png','-dpng');
% 
% figure(4);
% imagesc([0:7],[0:7],squeeze(stats.scorepredict(:,:,2)));
% aet;
% set(gca,'XAxisLocation','top');
% ylabel('NYCFC Goals');
% xlabel('Opponent Goals');
% title('Draw Probability');
% print('./output/DrawProbability.png','-dpng');
% 
% figure(5);
% imagesc([0:7],[0:7],squeeze(stats.scorepredict(:,:,3)));
% aet;
% set(gca,'XAxisLocation','top');
% ylabel('NYCFC Goals');
% xlabel('Opponent Goals');
% title('Loss Probability');
% print('./output/LossProbability.png','-dpng');


[gfs] = hist(stats.gfs,.001:5:95);
[gas] = hist(stats.gas,.001:5:95);
plot(gfs,'b');
hold on;
plot(gas,'r')
hold off

winind = 1:ngames;
winind = winind([game(:).win]);
[gfs] = hist(stats.gfs(ismember(stats.gfinds,winind)),.001:5:95);
[gas] = hist(stats.gas(ismember(stats.gainds,winind)),.001:5:95);
plot(gfs,'b');
hold on;
plot(gas,'r')
hold off

drawind = 1:ngames;
drawind = drawind([game(:).draw]);
[gfs] = hist(stats.gfs(ismember(stats.gfinds,drawind)),.001:5:95);
[gas] = hist(stats.gas(ismember(stats.gainds,drawind)),.001:5:95);
plot(gfs,'b');
hold on;
plot(gas,'r')
hold off

lossind = 1:ngames;
lossind = lossind([game(:).loss]);
[gfs] = hist(stats.gfs(ismember(stats.gfinds,lossind)),.001:5:95);
[gas] = hist(stats.gas(ismember(stats.gainds,lossind)),.001:5:95);
plot(gfs,'b');
hold on;
plot(gas,'r')
hold off

