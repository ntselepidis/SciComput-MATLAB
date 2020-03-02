function sep = vertex_separators(G,map,Out,couples)
ndoms=length(Out);
sep=cell(ndoms,1);
for i=1:size(couples,1)
    aborder=[];
    for j=Out{couples(i,1)}
        nbr=find(G(:,j));
        if any( map(nbr)==couples(i,2) )
            aborder(end+1)=j;
        end
    end
    bborder=[];
    for j=Out{couples(i,2)}
        nbr=find(G(:,j));
        if any( map(nbr)==couples(i,1) )
            bborder(end+1)=j;
        end
    end
    [p1,q1,r1,s1] = dmperm( G(aborder,bborder) );
    sepa = aborder( p1( r1(1):(r1(2)-1)) );
    sepb = bborder( q1( s1(2):(s1(length(s1))-1) ) );
    stemp = [sepa sepb];
    sep{couples(i,1)} = union( sep{couples(i,1)}, stemp );
    sep{couples(i,2)} = union( sep{couples(i,2)}, stemp );
end
for i=1:ndoms
    sep{i} = reshape( sep{i}, 1, numel(sep{i}) );
end
end