% Parameters

clear all;

n = 50;     % number of agents
d = 10;     % dimension

maxiter = 400; % maximal iteration number

%% objective functions

lambda = 0.001; %% coefficient of the regularization term
A = zeros(d*n,d*n);
b = zeros(d*n,1);

for ii = 1:n
    A((ii-1)*d+1:ii*d,(ii-1)*d+1:ii*d) = eye(d) + randn(d,d) * 0.1;
    b((ii-1)*d+1:ii*d,1) = randn(d,1); 
end

%% initial condition
x0 = zeros(d,n);

%% Generate graph
Adj = zeros(n,n); % adjacency matrix
pp = 0.1; % probability of generating an edge 


edges = {};
kkk = 10;
ll = 0;

for iii = 1:kkk
    for i = 1:n
        for  j = 1:n
            tmp = rand;
            if i > j && tmp >= 1-pp; 
                Adj(i,j) = 1;
                Adj(j,i) = 1;
                
                ll = ll + 1;
                edges{ll} = [j,i];
            end
        end
    end
    eig_lap = eig(laplacian(graph(Adj)));
    if eig_lap(2) >0 %% connectivity
        break
    end
end

%% generate the graph
G = graph(Adj);

mixing_matrix = eye(n) - 1/max( Adj * ones(n,1) ) * laplacian(G);

D = []; %% CD matrix for maximal cliques
DD = {}; %% {D_1, D_2, ....}

cliques = maximalCliques(Adj);

for l = 1:length(cliques)
    tmp = zeros(length(cliques{l}),n);

    for j = 1:length(cliques{l})
        tmp(j, cliques{l}(j) ) = 1; 
    end
    D = [D; tmp];
    DD{l} = tmp;
end

D_e = []; %% CD matrix for edges
DD_e = {}; %% {D_1, D_2, ... } for D_e

for l = 1:length(edges)
    tmp = zeros(length(edges{l}),n);

    for j = 1:length(edges{l})
        tmp(j, edges{l}(j) ) = 1; 
    end

    D_e = [D_e; tmp];
    DD_e{l} = tmp;
end

%% compute optimimal solution via CVX
cvx_begin
    variable xx(d)
    OBJ = 0;
    for ii = 1:n
        OBJ = OBJ + obj_quad(A((ii-1)*d+1:ii*d,(ii-1)*d+1:ii*d),b((ii-1)*d+1:ii*d,1),xx);
    end
    OBJ = OBJ + n* lambda * norm(xx,1);

    minimize OBJ
cvx_end

x_opt = xx;
f_opt = 0;

f_opt = obj_quad(A,b,kron(ones(n,1),x_opt));
f_opt = f_opt + n* lambda * norm(x_opt,1);


%% set an initial value
for i = 1:n
    x0(:,i) = randn(d,1)*1000;
    % x0(:,i) = randn(d,1)*1;
    % x0(:,i) = x_opt;
end