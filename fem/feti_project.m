function y = feti_project(x,FETI)
    if isempty(FETI.G)
        y = x;
        return;
    end
    G = FETI.G;
    y = x - G * ((G'*G)\(G'*x));
end

