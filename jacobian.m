function g = jacobian(r,x)
% g = jacobian(r,x)
%
% Calculates the jacobian of the vector valued function r at x.

lx = length(x);
g = [];
h = eps^(1/3);
h2=h*2;
for i = 1:lx
   xplus = x;
   xminus = x;
   xplus(i) = x(i) + h;
   xminus(i) = x(i) - h;
   g(:,i) = ( r(xplus) - r(xminus) )/h2;
end
