function overlap = form_overlap_ml(G, In, Out, delta)

    n       = length(G);
    ndoms   = length(Out);
    overlap = cell(ndoms, 1);

    for i = 1 : ndoms
        v = zeros(n, 1);
        v( Out{i} ) = 1;
        for lev = 1 : delta
            ovpnts = find( v == lev );
            for j = 1 : length(ovpnts)
                pnt = ovpnts(j);
                nbr = find( G(:,pnt) );
                v(nbr) = lev + 1;
            end
            if ( lev == 1 )
                v( In{i} ) = 0;
            end
        end
        overlap{i} = find( v > 0 );
        overlap{i} = reshape( overlap{i}, 1, length(overlap{i}) );
    end
    
end
