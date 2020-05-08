function [nbr_doms, T, couples] = neighbor_domains(G, map, Out)
ndoms = length(Out);
nbr_doms = cell(ndoms,1);
for i = 1 : ndoms
    temp_doms = zeros(1,ndoms);
    for j = Out{i}
        nbr = find(G(:,j));
        map_neq_i = ( map(nbr) ~= i );
        nbr_doms_j = map( nbr(map_neq_i) );
%         nbr_doms{i} = union(nbr_doms{i}, nbr_doms_j);
        temp_doms(nbr_doms_j) = 1;
    end
    nbr_doms{i} = find(temp_doms);
end
Ti = [];
Tj = [];
for i = 1 : ndoms
    for j = nbr_doms{i}
        Ti(end+1) = i;
        Tj(end+1) = j;
    end
end
T = sparse(Ti, Tj, ones(1,length(Ti)), ndoms, ndoms);
T = tril(T);
couples = zeros(nnz(T),2);
idx = 1;
for i = 1 : ndoms
    nbr = find(T(:,i));
    for j = 1 : length(nbr)
        couples(idx,:) = [i nbr(j)];
        idx = idx+1;
    end
end
end
