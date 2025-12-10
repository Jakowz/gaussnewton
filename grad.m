function g = grad(f,x)
% g = grad(f,x)
%
% Calculates the gradient (column) vector of the function f at x.

lx = length(x);
g = zeros(lx,1);
h = eps^(1/3);
h2=h*2;
for i = 1:lx
   xplus = x;
   xminus = x;
   xplus(i) = x(i) + h;
   xminus(i) = x(i) - h;
   g(i,1) = ( f(xplus) - f(xminus) )/h2;
end
