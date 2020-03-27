function [X] = ConstructTensorTri(dm,mm)

    [nd,nm]=size(dm);
    X=zeros(nm,nm,nd);
    for n=1:nd
        i=find(dm(n,:)==1);
        for m=1:length(i)
            j=find(mm(i(m),:)==1);        
            s=intersect(j,i);
            X(i(m),s,n)=1;
        end

    end

    X=tensor(X);
end

