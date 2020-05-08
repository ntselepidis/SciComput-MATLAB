function [R, basis] = shem(A, In, Out, All, nbasis)

    n      = length(A);
    ndoms  = length(In);
    basis  = cell(ndoms, 1);
    source = cell(ndoms, 1);
    Sc     = cell(ndoms, 1);
    for i=1:ndoms

        ninn = length(In{i} );
        nout = length(Out{i});
        basis{i} = zeros(ninn+nout, nbasis);    
        basis{i}(:,1) = 1;

        if ( nbasis > 1 )

            Sc{i} = A(Out{i},Out{i}) - A(Out{i},In{i}) * ( A(In{i},In{i}) \ A(In{i},Out{i}) );
%             [source{i}, ~] = eig( full( Sc{i} ) );
%             [source{i}, ~] = eigs( Sc{i}, nbasis-1, 'smallestabs' );
            [source{i}, ~] = eigs( Sc{i}, nbasis-1, 'SM' );

            for j = 2 : nbasis 
                basis{i}(ninn+1:end,j) = source{i}(:,j-1);
                basis{i}(1:ninn,    j) = A(In{i},In{i}) \ ( - A(In{i},Out{i}) * basis{i}(ninn+1:end,j) );
            end

        end

    end
    
    R = spalloc(nbasis*ndoms, n, nbasis*n);
    for i = 1 : ndoms
        for j = 1 : nbasis
            R(nbasis*(i-1)+j, All{i}) = basis{i}(:,j)';
        end
    end

end
