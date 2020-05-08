function [In, Out, All, blk] = inner_outer(G, map)

    ndoms = max(map);
    In  = cell(ndoms, 1);
    Out = cell(ndoms, 1);
    All = cell(ndoms, 1);
    blk = zeros(ndoms+1, 1);
    blk(1) = 1;
    for i = 1 : ndoms
        dpnts = find( map == i );
        for j = 1 : length(dpnts)
            pnt = dpnts(j);
            nbr = find(G(:,pnt));
            if any( map(nbr) ~= i )
                Out{i} = [Out{i} pnt];
            else
                In{i} = [In{i} pnt];
            end
        end
        All{i} = [In{i} Out{i}];
        blk(i+1) = blk(i) + length(All{i});
    end

end
