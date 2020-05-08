function [s, deg, p] = global_separators(All, sep)
p = horzcat(All{:});
s2 = horzcat(sep{:});
ss = sort(s2);
s = unique(ss);
deg = ones(1,length(s));
idx = 1;
t1 = ss(1);
for i = 2 : length(ss)
    if ( ss(i) == t1 )
        deg(idx) = deg(idx)+1;
    else
        t1 = ss(i);
        idx = idx+1;
    end
end
end
