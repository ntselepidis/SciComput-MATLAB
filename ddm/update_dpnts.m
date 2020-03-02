function [In, All, blk] = update_dpnts(All, sep, n)
ndoms=length(sep);
In=cell(ndoms,1);
blk=zeros(ndoms+1,1);
blk(1)=1;
for i=1:ndoms
    ind=zeros(1,n);
    ind( All{i} ) = 1; 
    ind( sep{i} ) = 0;
    In{i}=find(ind);
    All{i}=[In{i} sep{i}];
%     All{i}=[setdiff(All{i},sep{i}) sep{i}];
    blk(i+1)=blk(i)+length(All{i});
end
end