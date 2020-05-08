function [All2, overlap] = form_overlap(G, map, All, Out)
    
    n       = length(G);
    ndoms   = length(All);
    overlap = cell(ndoms, 1);
    All2    = cell(ndoms, 1); % All including overlap
    
    % one level overlap (delta = 1)
    for i = 1 : ndoms
        pnts = zeros(n, 1);
        for j = Out{i}
            nbr = find(G(:,j));
            map_neq_i = ( map(nbr) ~= i );
            pnts( nbr(map_neq_i) ) = 1;
        end
        overlap{i} = find(pnts)';
        All2{i} = [All{i} overlap{i}];
    end
    
end
