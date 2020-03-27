function [X_new,C_tensor,U,S,V,errX] = miRCom(O, W, A, B,LA,LB,Ui,Vj,maxIter,K,R,alpha,beta)

n_c = length(A);     
n_d = length(B);      
tol = 1e-4;           

% Lagrange multipliers
rho = 1e-6;
eta = 1e-6;
mu = 1e-6;
gama= 1e-6;
epsilon = 1e-6;

rhomax = 1e6;
etamax = 1e6;
mumax = 1e6;
gamamax = 1e6;
epsilonmax = 1e6;

coeff = 1.15; 
errX = zeros(maxIter, 1);

% compute statistics for O
dim = size(O); 

%% initialization

X = O; 
U = randn(dim(1),K);
S = U;
V = randn(dim(3),K); 

% side info drugs
Q = cell(1,n_c);
for i = 1:n_c
    Q{i} = randn(K,R);
end
C = Ui;
H = Ui;

% side info  diseases
P = cell(1,n_d);
for j = 1:n_d
    P{j}=randn(K,R);
end
D = Vj;
G = Vj;


% Lagrange multipliers
F = zeros(dim(1), K);

Y = cell(1,n_c);
E=cell(1,n_c);
for i = 1:n_c
    Y{i} = zeros(dim(1),R);
    E{i} = zeros(dim(1),R);
end

Z = cell(1,n_d);
L=cell(1,n_d);
for j = 1:n_d
    Z{j} = zeros(dim(3),R);
    L{j} = zeros(dim(3),R);
end


%% main loop

for iter = 1: maxIter       
    if mod(iter, 20) == 0
        fprintf('miRCom: iterations = %d   difference=%f\n', iter, errX(iter-1));
    end

    % update the U
    Pi = 2* (V' * V) .*(S' * S) + rho * eye(K) ;    
    c = 2 * mttkrp(X,{U,S,V},1) + rho * S + F;

    for i = 1:n_c
        c = c + 2 * alpha * Ui{i} * Q{i}';
        Pi=Pi+ 2 * alpha *  Q{i}* Q{i}';
    end
    U = c / Pi;

    %Update the S
    Theta = 2 * (V' * V) .*(U' * U) + rho * eye(K);
    c = 2 * mttkrp(X,{U,S,V},2) + rho * U - F;
    S = c / Theta;


    % Update the V
    Xi = (S' * S) .*(U' * U);
    c = mttkrp(X,{U,S,V},3);
    
    for j = 1:n_d
        c = c + beta * Vj{j} * P{j}';
        %temp = size(LB{j},1);
        Xi = Xi + beta * P{j} * P{j}';
    end
    V = c / Xi;


    % Update the Ui n*R
    for i = 1:n_c
        a = 2 * A{i}' * C{i} + 2 * U * Q{i} + eta * C{i} + gama * H{i} + Y{i}+E{i};
        b = 2 * C{i}' * C{i} + 2 * eye(R) + eta * eye(R)+ gama * eye(R);
        Ui{i} = a / b;
        %Q{i} = diag(sum(Ui{i}));
    end 


    % Update the C{i} n*R
    for i = 1:n_c
        a = 2 * A{i} * Ui{i} + eta * Ui{i} - Y{i};
        b = 2 * Ui{i}' * Ui{i}  + eta * eye(R);
        C{i} = a / b; 
    end 

    
    %Update the Q{i}
    for i=1:n_c
       temp=Q{i};
       len=size(temp,1);
       temp2=zeros(len,len);
       for j=1:len
            temp2(j,j)=1/2*norm(temp(j,:),2);
       end
       tri=temp2;
       a= U' * Ui{i};
        b=(U' * U) + tri;
        Q{i}=b^-1 * a;
    end
    
    %Update the H{i} n*R
    for i = 1:n_c
        a = gama * Ui{i} - E{i};
        temp=size(LA{i},1);
        b = 2 * LA{i}  + gama * eye(temp);
        H{i} = b^-1*a; 
    end 

   % Update the Vj
    for j = 1:n_d
        a = 2 * B{j}' * D{j} + 2 * V * P{j} + mu * D{j} + Z{j} + epsilon * G{j} + L{j};
        b = 2 * D{j}' * D{j} + 2 * eye(R) + mu * eye(R) + epsilon *eye(R);
        Vj{j} = a / b;
        %P{j} = diag(sum(Vj{j}));
    end 


    % Update the D{j}
    for  j= 1:n_d
        a = 2 * B{j} * Vj{j} + mu * Vj{j} - Z{j};
        b = 2 * Vj{j}' * Vj{j}  + mu * eye(R);
        D{j} = a / b;
    end 
    
    %Update the P{i}
    for i=1:n_d
       temp=P{i};
       len=size(temp,1);
       temp2=zeros(len,len);
       for j=1:len
            temp2(j,j)=1/2*norm(temp(j,:),2);
       end
       tri=temp2;
       a= V' * Vj{i};
        b=(V' * V) + tri;
        P{i}=b^-1 * a;
    end
    
    %Update the G{j} m*R
    for i = 1:n_d
        a = epsilon * Vj{i} - L{i};
        temp=size(LB{i},1);
        b = 2 * LB{i}  + epsilon * eye(temp);
        G{i} = (b^-1)*a; 
    end 
    
    % Update the Lagrange Multipler
    F = F + rho *(S - U);

    for i = 1:n_c
        Y{i} = Y{i} + eta * (C{i} - Ui{i});
        E{i} = E{i} + gama * (H{i} - Ui{i});
    end

    for j = 1:n_d
        Z{j} = Z{j} + mu * (D{j} - Vj{j});
        L{j} = L{j} + epsilon * (G{j} - Vj{j});
    end
    

    % update mu, eta, rho
    rho = min(coeff * rho, rhomax);
    mu = min(coeff * mu, mumax);
    eta = min(coeff * eta, etamax);
    gama = min(coeff * gama, gamamax);
    epsilon = min(coeff * epsilon, epsilonmax);

    % compute the fit.
    C_tensor = ktensor({U,S,V});
%   C_tensor=tensor(C_tensor);
%   X_new=double(C_tensor);
%   X_new(W) = O(W);

    % update X
    TEMP = times(1-W,C_tensor); 
    X_new = plus(full(O),full(TEMP));
    errX(iter) = norm(X(:) - X_new(:))/norm(X(:)) ; 

    X=X_new; 

    % Check for convergence
    if (iter > 1) && (errX(iter) < tol)
        break;
    end
   
end
    C_tensor=tensor(C_tensor);
end

