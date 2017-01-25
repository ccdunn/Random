%albania
ldw = [0 1 3];

A = zeros(3,3,3,3,3,3);
B = zeros(3,3,3,3,3,3);
C = zeros(3,3,3,3,3,3);
D = zeros(3,3,3,3,3,3);
ii = 0;
for aa = 1:3 %A vs B
    for bb = 1:3 %A vs C
        for cc = 1:3 %A vs D
%             A(aa,bb,cc,:,:,:) = ldw(aa) + ldw(bb) + ldw(cc);
            for dd = 1:3 %B vs C
                for ee = 1:3 %B vs D
%                     if bb==1 && cc==1
%                         B(aa,:,:,dd,ee,:) = ldw(4 - aa) + ldw(dd) + ldw(ee);
%                     end

                    for ff = 1:3 %C vs D
%                         if aa==1 && cc==1 && ee==1
%                             C(:,bb,:,dd,:,ff) = ldw(4 - bb) + ldw(4 - dd) + ldw(ff);
%                         end

%                         if aa==1 && bb==1 && dd==1
%                             D(:,:,cc,:,ee,ff) = ldw(4 - cc) + ldw(4 - ee) + ldw(4 - ff);
%                         end
                        A(aa,bb,cc,dd,ee,ff) = ldw(aa) + ldw(bb) + ldw(cc);
                            C(aa,bb,cc,dd,ee,ff) = ldw(4 - bb) + ldw(4 - dd) + ldw(ff);
                        B(aa,bb,cc,dd,ee,ff) = ldw(4 - aa) + ldw(dd) + ldw(ee);
                            D(aa,bb,cc,dd,ee,ff) = ldw(4 - cc) + ldw(4 - ee) + ldw(4 - ff);
                    end
                end
            end
        end
    end
end

A = A(:);
B = B(:);
C = C(:);
D = D(:);

scores = cat(2,A,B,C,D);

[scores_s,ind] = sort(scores,2,'descend');

b = hist(scores_s(:,3),0:9);

b = b./sum(b);

stem(0:9,b);

for ii = 0:9
    pt = (sum(b(ii+2:end)) + b(ii+1)/2);
    p(ii+1) = 1 - pt.^4*(1-pt) - pt.^5;
end

stem(0:9,p);

probqualas3 = sum(b.*p);