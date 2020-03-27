function [normL,eigvecNormL] = MyLapalacian(sim,k)
    nd=size(sim,2);
    normL=cell(1,nd);
    eigvecNormL=cell(1,nd);
    for i=1:nd
        W=sim{i};
         s = sum(W,2); 
         %%normL
          s1=power(s,-0.5);
          D1=diag(s1);
          E=D1*W*D1;
          n=size(E,1);
          normL{i} = eye(n)-E;  
         [eigvecNormL{i},~] = MygetEigvec(normL{i},k);

    end

end

