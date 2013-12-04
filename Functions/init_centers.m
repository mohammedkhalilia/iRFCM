function V = init_centers(type, n, c, D)
    
    switch type
        case 1
            V = rand(c,n);
            V = V./(sum(V,2) * ones(1,n));
            
        case 2
            idx = randperm(n,c);
            V = D(idx,:);
            V = V./(sum(V,2) * ones(1,n));
    end
end