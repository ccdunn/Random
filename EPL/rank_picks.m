
clear all;
fn = fullfile(homedir,'personal','EPL','EPL_15-16_entries.csv');
tfn = fullfile(homedir,'personal','EPL','EPL_15-16_truth.csv');

addpath(genpath(fullfile(homedir,'libs','stringdist')));
addpath(genpath(fullfile(homedir,'libs','munkres')));

fid = fopen(fn);
entries_in_txt = textscan(fid,'%s','delimiter','\n\r');
fclose(fid);
entries_in_txt = entries_in_txt{1};

fid = fopen(tfn);
truth_in_txt = textscan(fid,'%s','delimiter','\n\r');
fclose(fid);
truth_in_txt = truth_in_txt{1};

entries_txt = {};
for ii = 1:size(entries_in_txt,1)
    entries_txt(ii,:) = strsplit(entries_in_txt{ii},',');
end
truth_txt = {};
for ii = 1:size(entries_in_txt,1)
    truth_txt(ii,:) = strsplit(truth_in_txt{ii},',');
end

names = entries_txt(1,:);
entries_in = entries_txt(2:end,:);
truth_in = truth_txt(2:end,:);



P = numel(names);
N = size(entries_in,1);

lookup = {...
    {'AFC Bournemouth','Bournemouth'},...
    {'Arsenal'},...
    {'Aston Villa','Villa'},...
    {'Chelsea'},...
    {'Crystal Palace'},...
    {'Everton'},...
    {'Leicester City','Leicester'},...
    {'Liverpool'},...
    {'Manchester City','ManCity'}...
    {'Manchester United','ManU'},...
    {'Newcastle United','Newcastle'},...
    {'Norwich City','Norwich'},...
    {'Southampton'},...
    {'Stoke City','Stoke'},...
    {'Sunderland'},...
    {'Swansea City','Swansea'}...
    {'Tottenham Hotspur','Spurs'},...
    {'Watford'},...
    {'West Bromwich Albion','WBA'},...
    {'West Ham United','West Ham'}};

lookup_p = lookup;
for ii = 1:N
    lookup_p{ii} = lookup_p{ii}{1};
end

w = zeros(N,N);
entries = zeros(size(entries_in));
for pp = 1:P
    
    for ii = 1:N
        for jj = 1:N
            m = Inf(1,length(lookup{ii}));
            for kk = 1:length(lookup{ii})
                m(kk) = strdist(lookup{ii}{kk},entries_in{jj,pp});
            end
            w(jj,ii) = min(m);
        end
    end
    
    entries(:,pp) = munkres(w);
    
end


truth = truth_in;
truth = zeros(size(truth));
for ii = 1:N
    for jj = 1:N
        m = Inf(1,length(lookup{ii}));
        for kk = 1:length(lookup{ii})
            m(kk) = strdist(lookup{ii}{kk},truth_in{jj});
        end
        w(jj,ii) = min(m);
    end
end

truth = colvec(munkres(w));

%reorganize so place listed by team, not team listed by place
ind = colvec(1:N);

truth0 = truth;
truth(truth) = ind;
entries0 = entries;
for pp = 1:P
    entries(entries(:,pp),pp) = ind;
end

truth
entries

scores = zeros(N,P);
score = zeros(P,1);
for pp = 1:P
    [score(pp),scores(:,pp)] = scorefun(truth,entries(:,pp));
end

[score,ranks] = sort(score,'ascend');

scores = -scores(:,ranks);
scores = scores(truth0,:);

minv = repmat(min(abs(scores),[],2),1,P);
mins = minv==abs(scores);
[minrow,mincol] = find(mins);

labs = cell(1,P+1);
for pp = 1:P
    labs{pp} = sprintf('%s (%i)',names{ranks(pp)},score(pp));
end
labs{P+1} = 'Mean';

c1 = .05;
c2 = .8;

cmap1 = cat(2,colvec(linspace(c1,c2,N)),colvec(linspace(c1,c2,N)),ones(N,1));
cmap2 = cat(2,ones(N,1),colvec(linspace(c2,c1,N)),colvec(linspace(c2,c1,N)));
cmap = cat(1,cmap1,[1 1 1],cmap2);

m = mean(scores,2);

imagesc(cat(2,scores,m));
hold on;
scatter(mincol,minrow,'k.');
hold off;
aet;
caxis([-20 20]);
cbar = colorbar;
colormap(cmap);
set(gca,'YTick',1:N);
set(gca,'YTickLabel',lookup_p(truth0));
set(gca, 'XAxisLocation', 'top')
set(gca,'XTick',1:P+1);
set(gca,'XTickLabel',labs);
set(gca,'XTickLabelRotation',90);
set(cbar,'YTick',[]);
text(P+3,10.5 - 3,'Underestimated','Rotation',90,'HorizontalAlignment','left');
text(P+3,10.5 + 3,'Overestimated','Rotation',90,'HorizontalAlignment','right');
title('Results (L1 Norm)');



