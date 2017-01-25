r = 1:9;

[a,b] = meshgrid(r,r);

no = sort([colvec(a), colvec(b)],2,'ascend');
no = unique(no,'rows');
oind = (1:length(no)).';

pat = zeros(length(no),1);
len = {};

count = 1;

po = prod(no,2);
so = sum(no,2);

n = no;
p = po;
s = so;
nind = oind;

nindprev = [-1];

while(~isempty(nind))
    
    if mod(count,2)
        x = p;
    else
        x = s;
    end
    
    b = min(x(:)):max(x(:));
    nx = hist(x,b);
    
    indn = ismember(x,b(nx>1));
    indnd = ismember(x,b(nx==1));
    
    len{count} = oind(nind(indnd));
    pat(oind(nind(indnd))) = count;
    p = p(indn);
    s = s(indn);
    n = n(indn,:);
    nind = nind(indn,:);   
    
    count = count + 1;
    
    if numel(nind)==numel(nindprev) && all(nind==nindprev)
        break;
    end 
    nindprev = nind;
end