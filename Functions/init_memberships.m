function U = init_memberships(type, n, c)
    
    switch type
        case 1
            psize=ceil(n/c);
            temp=1;
            for i=1:n
                ind(i)=temp;
                if mod(i,psize)==0
                    temp=temp+1;
                end
            end

            %produce the c x n partition matrix
            U = sparse(ind,1:n,1,c,n);
        case 2
            U = rand(c,n);
            U = U./(ones(c,1)*sum(U));
    end
end