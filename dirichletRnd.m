function x = dirichletRnd(n, m)
if (nargin == 2)
    n = n*m;
end
x = gamrnd(n,1);
x = x / sum(x);
end