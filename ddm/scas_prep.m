function AS = scas_prep(A, SC, overlap, s_droptol, ilu_setup)

if nargin == 4
    inc = 0;
    disp('complete schurAS');
else
    inc = 1;
    disp('incomplete schurAS');
end

% ATTENTION : ONLY ONE-LEVEL OVERLAP IS CURRENTLY SUPPORTED
[In, Out, L, U, P, Q] = deal(SC.In, SC.Out, SC.L, SC.U, SC.P, SC.Q);
ndoms = length(In);
[S, Ls, Us, Ps, Qs, T] = deal( cell(ndoms,1) );
for i=1:ndoms
    overlap{i} = [Out{i} setdiff(overlap{i},Out{i})]; % REORDER OVERLAP VERTICES

    if (inc == 0)
        T{i} = A(overlap{i},In{i}) * sp_solve( L{i},U{i},P{i},Q{i}, A(In{i},overlap{i}) );
    else
        [iL,iU] = ilu( A(In{i},In{i}), ilu_setup );
        T{i} = A(overlap{i},In{i}) * ( iU \ ( iL \ (A(In{i},overlap{i})) ) );
    end
    
    S{i} = A(overlap{i},overlap{i}) - T{i};
    
    % sparsify S
    
    D = abs(diag(S{i}));
    
    if ( s_droptol ~= 0 )
        for r = 1 : size(S{i},1)
            for c = 1 : size(S{i},2)
                if ( abs(S{i}(r,c)) <= s_droptol * (D(r)+D(c)) )
                    S{i}(r,c) = 0;
                end
            end
        end
    end
    
    [Ls{i}, Us{i}, Ps{i}, Qs{i}] = lu( S{i}, 'vector' );
end
[AS.L, AS.U, AS.P, AS.Q] = deal( Ls, Us, Ps, Qs );
AS.T = T;

% MASKING OVERLAP BASED ON OUTER POINTS (ONLY FOR DELTA = 1)
out = horzcat(Out{:});
mask( out ) = 1 : length(out);
for i = 1 : ndoms
    overlap{i} = mask( overlap{i} );
end
AS.overlap = overlap;

end
