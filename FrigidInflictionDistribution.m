
% items to distribute
items = {'mug','socks','mammut hat','tshirt',...
    'knife','headlamp','deodorizer','OR hat',...
    'backpack','cider green hat','tick twister','magazine',...
    'chalk bag'};

% everyone's ranking of the objects. Subtracted from 13 so that more
% desired objects are a higher value
J = 14 - [ 1  2  3  4  5  6  7  8  9 10 11 12 13];
A = 14 - [ 4  6  1  2  3  9 11  7  5 12 10 13  8];
C = 14 - [ 5  1  9 12  4  6  3  7  2 13  8 11 10];

J = 2.^(J);
A = 2.^(A);
C = 2.^(C);


%loop through all possible combinations (only works for ~13 items and ~3 people)
J_score = zeros(3,3,3,3,3,3,3,3,3,3,3,3,3);
A_score = zeros(3,3,3,3,3,3,3,3,3,3,3,3,3);
C_score = zeros(3,3,3,3,3,3,3,3,3,3,3,3,3);

for aa = 1:3
    J_score(aa,:,:,:,:,:) = J_score(aa,:,:,:,:,:) + J(1)*(aa==1);
    A_score(aa,:,:,:,:,:) = A_score(aa,:,:,:,:,:) + A(1)*(aa==2);
    C_score(aa,:,:,:,:,:) = C_score(aa,:,:,:,:,:) + C(1)*(aa==3);
    for bb = 1:3
        J_score(aa,bb,:,:,:,:) = J_score(aa,bb,:,:,:,:) + J(2)*(bb==1);
        A_score(aa,bb,:,:,:,:) = A_score(aa,bb,:,:,:,:) + A(2)*(bb==2);
        C_score(aa,bb,:,:,:,:) = C_score(aa,bb,:,:,:,:) + C(2)*(bb==3);
        for cc = 1:3
            J_score(aa,bb,cc,:,:,:) = J_score(aa,bb,cc,:,:,:) + J(3)*(cc==1);
            A_score(aa,bb,cc,:,:,:) = A_score(aa,bb,cc,:,:,:) + A(3)*(cc==2);
            C_score(aa,bb,cc,:,:,:) = C_score(aa,bb,cc,:,:,:) + C(3)*(cc==3);
            for dd = 1:3
                J_score(aa,bb,cc,dd,:,:) = J_score(aa,bb,cc,dd,:,:) + J(4)*(dd==1);
                A_score(aa,bb,cc,dd,:,:) = A_score(aa,bb,cc,dd,:,:) + A(4)*(dd==2);
                C_score(aa,bb,cc,dd,:,:) = C_score(aa,bb,cc,dd,:,:) + C(4)*(dd==3);
                for ee = 1:3
                    J_score(aa,bb,cc,dd,ee,:) = J_score(aa,bb,cc,dd,ee,:) + J(5)*(ee==1);
                    A_score(aa,bb,cc,dd,ee,:) = A_score(aa,bb,cc,dd,ee,:) + A(5)*(ee==2);
                    C_score(aa,bb,cc,dd,ee,:) = C_score(aa,bb,cc,dd,ee,:) + C(5)*(ee==3);
                    for ff = 1:3
                        J_score(aa,bb,cc,dd,ee,ff) = J_score(aa,bb,cc,dd,ee,ff) + J(6)*(ff==1);
                        A_score(aa,bb,cc,dd,ee,ff) = A_score(aa,bb,cc,dd,ee,ff) + A(6)*(ff==2);
                        C_score(aa,bb,cc,dd,ee,ff) = C_score(aa,bb,cc,dd,ee,ff) + C(6)*(ff==3);
                        for gg = 1:3
                            J_score(aa,bb,cc,dd,ee,ff,gg) = J_score(aa,bb,cc,dd,ee,ff,gg) + J(7)*(gg==1);
                            A_score(aa,bb,cc,dd,ee,ff,gg) = A_score(aa,bb,cc,dd,ee,ff,gg) + A(7)*(gg==2);
                            C_score(aa,bb,cc,dd,ee,ff,gg) = C_score(aa,bb,cc,dd,ee,ff,gg) + C(7)*(gg==3);
                            
                            for hh = 1:3
                                J_score(aa,bb,cc,dd,ee,ff,gg,hh) = J_score(aa,bb,cc,dd,ee,ff,gg,hh) + J(8)*(hh==1);
                                A_score(aa,bb,cc,dd,ee,ff,gg,hh) = A_score(aa,bb,cc,dd,ee,ff,gg,hh) + A(8)*(hh==2);
                                C_score(aa,bb,cc,dd,ee,ff,gg,hh) = C_score(aa,bb,cc,dd,ee,ff,gg,hh) + C(8)*(hh==3);
                                
                                for ii = 1:3
                                    J_score(aa,bb,cc,dd,ee,ff,gg,hh,ii) = J_score(aa,bb,cc,dd,ee,ff,gg,hh,ii) + J(9)*(ii==1);
                                    A_score(aa,bb,cc,dd,ee,ff,gg,hh,ii) = A_score(aa,bb,cc,dd,ee,ff,gg,hh,ii) + A(9)*(ii==2);
                                    C_score(aa,bb,cc,dd,ee,ff,gg,hh,ii) = C_score(aa,bb,cc,dd,ee,ff,gg,hh,ii) + C(9)*(ii==3);
                                    
                                    for jj = 1:3
                                        J_score(aa,bb,cc,dd,ee,ff,gg,hh,ii,jj) = J_score(aa,bb,cc,dd,ee,ff,gg,hh,ii,jj) + J(10)*(jj==1);
                                        A_score(aa,bb,cc,dd,ee,ff,gg,hh,ii,jj) = A_score(aa,bb,cc,dd,ee,ff,gg,hh,ii,jj) + A(10)*(jj==2);
                                        C_score(aa,bb,cc,dd,ee,ff,gg,hh,ii,jj) = C_score(aa,bb,cc,dd,ee,ff,gg,hh,ii,jj) + C(10)*(jj==3);
                                        
                                        for kk = 1:3
                                            J_score(aa,bb,cc,dd,ee,ff,gg,hh,ii,jj,kk) = J_score(aa,bb,cc,dd,ee,ff,gg,hh,ii,jj,kk) + J(11)*(kk==1);
                                            A_score(aa,bb,cc,dd,ee,ff,gg,hh,ii,jj,kk) = A_score(aa,bb,cc,dd,ee,ff,gg,hh,ii,jj,kk) + A(11)*(kk==2);
                                            C_score(aa,bb,cc,dd,ee,ff,gg,hh,ii,jj,kk) = C_score(aa,bb,cc,dd,ee,ff,gg,hh,ii,jj,kk) + C(11)*(kk==3);
                                            
                                            for ll = 1:3
                                                J_score(aa,bb,cc,dd,ee,ff,gg,hh,ii,jj,kk,ll) = J_score(aa,bb,cc,dd,ee,ff,gg,hh,ii,jj,kk,ll) + J(12)*(ll==1);
                                                A_score(aa,bb,cc,dd,ee,ff,gg,hh,ii,jj,kk,ll) = A_score(aa,bb,cc,dd,ee,ff,gg,hh,ii,jj,kk,ll) + A(12)*(ll==2);
                                                C_score(aa,bb,cc,dd,ee,ff,gg,hh,ii,jj,kk,ll) = C_score(aa,bb,cc,dd,ee,ff,gg,hh,ii,jj,kk,ll) + C(12)*(ll==3);
                                                
                                                for mm = 1:3
                                                    J_score(aa,bb,cc,dd,ee,ff,gg,hh,ii,jj,kk,ll,mm) = J_score(aa,bb,cc,dd,ee,ff,gg,hh,ii,jj,kk,ll,mm) + J(13)*(mm==1);
                                                    A_score(aa,bb,cc,dd,ee,ff,gg,hh,ii,jj,kk,ll,mm) = A_score(aa,bb,cc,dd,ee,ff,gg,hh,ii,jj,kk,ll,mm) + A(13)*(mm==2);
                                                    C_score(aa,bb,cc,dd,ee,ff,gg,hh,ii,jj,kk,ll,mm) = C_score(aa,bb,cc,dd,ee,ff,gg,hh,ii,jj,kk,ll,mm) + C(13)*(mm==3);
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

J_score = J_score./sum(J);
A_score = A_score./sum(A);
C_score = C_score./sum(C);

step = .001;
lambdas = step:step:1 - step;
for ll = 1:length(lambdas)
    %find mean score from each combination
    mean_score = (J_score + A_score + C_score)/3;
    
    %find variance from each combination
    var_score = sqrt((J_score - mean_score).^2 + (A_score - mean_score).^2 + (C_score - mean_score).^2);
    
    %now, we want a solution that is a combination of high mean score and small
    %variance (so most overall happiness with least envy)
    %use lambda to weight between these two different criteria
    %high lambda weights overal happiness, low lambda weights fairness
    lambda = lambdas(ll);
    %lambda = .00000001 and .9999999 are interesting
    fs = -lambda*mean_score + (1-lambda)*var_score;
    
    %find minimum of cost function
    [m,mi] = min(fs(:));
    [fa(1),fa(2),fa(3),fa(4),fa(5),fa(6),fa(7),fa(8),fa(9),fa(10),fa(11),fa(12),fa(13)] = ind2sub(size(J_score),mi);
    
    J_items{ll} = items(find(fa==1));
    A_items{ll} = items(find(fa==2));
    C_items{ll} = items(find(fa==3));
    
    Js(ll) = J_score(mi);
    As(ll) = A_score(mi);
    Cs(ll) = C_score(mi);
    ms(ll) = mean_score(mi);
    vs(ll) = var_score(mi);
end

plot(lambdas,Js,'r');
hold on;
plot(lambdas,As,'g');
plot(lambdas,Cs,'b');
plot(lambdas,ms,'k--');
plot(lambdas,vs,'k-.');
hold off;

xlabel('lambda');
ylabel('happiness');

legend({'Jeff','Ashin','Charlie','Mean','Variance'},'Location','Northwest');

