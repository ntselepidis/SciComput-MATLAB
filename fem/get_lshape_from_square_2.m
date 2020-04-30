function [coo, con, bounds, sep, dpnts] = get_lshape_from_square_2(a, b, m, coo1, con1)

n = (m+1)^2;
coo2 = zeros(size(coo1)); con2 = zeros(size(con1)); nodes2 = zeros(1,n);
coo3 = zeros(size(coo1)); con3 = zeros(size(con1)); nodes3 = zeros(1,n);

west = 1:(m+1):(n-m);
nodes2(west) = west+m; % west2 same as east1
for i = 1:m+1
    idx = (i-1)*(m+1)+2:i*(m+1); % attention to +2
    nodes2(idx) = idx-1+n-(i-1); % continue numbering
end

south = 1:(m+1);
nodes3(south) = south+(m+1)*m; % south3 same as north1
nodes3(m+2:n) = (2*n-m):(3*n-2*m-2); % continue numbering

idx = 1;
if (size(con1,2) == 4)
    for i = 1:m
        for j = 1:m
            con2(idx,:) = nodes2(con1(idx,:));
            con3(idx,:) = nodes3(con1(idx,:));
            idx = idx+1;
        end
    end
else % if (size(con1,2) == 3)
    for i = 1:m
        for j = 1:m
            con2(idx,:) = nodes2(con1(idx,:));
            con2(idx+1,:) = nodes2(con1(idx+1,:));
            con3(idx,:) = nodes3(con1(idx,:));
            con3(idx+1,:) = nodes3(con1(idx+1,:));
            idx = idx+2;
        end
    end
end
coo2(:,1) = coo1(:,1)+(b-a);
coo2(:,2) = coo1(:,2);
coo2(west,:) = [];

coo3(:,1) = coo1(:,1);
coo3(:,2) = coo1(:,2)+(b-a);
coo3(south,:) = [];

coo = [coo1;coo2;coo3];
con = [con1;con2;con3];

lsouth = [1:m+1 nodes2(2:m+1)];
lnorth = [nodes3((1:m+1)+(m+1)*m) nodes2((2:m+1)+(m+1)*m)];
lwest = [(m+2):(m+1):(n-2*m-1) nodes3(1:(m+1):(n-2*m-1))];
least = [nodes2(((m+2):(m+1):(n-2*m-1))+m) nodes3((1:(m+1):(n-2*m-1))+m)];
bounds = [lsouth lnorth lwest least];

sep = [nodes2(west) nodes3(south(1:end-1))];

dpnts = cell(3,1);
dpnts{1} = 1:n;
dpnts{2} = nodes2;
dpnts{3} = nodes3;

end
