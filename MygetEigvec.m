function [eigvec,eigval] = MygetEigvec(A, k)
% compute the k-leading eigenvector of A
    [v,d] = eig(A);
    d = diag(d);
    [~,idx] = sort(d);
    idx1 = idx(1:k);
    eigval = d(idx1);
    eigvec = v(:,idx1); 
 end