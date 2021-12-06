% Gradient of function f

function Gf = gradf(x,u)

Gf = zeros(4,4);
Gf(2,1) = 1;
Gf(1,2) = -4*x(1)*x(2)-1+x(2)*cos(x(1)*x(2))*(u+1/3*u^3+sin(u));
Gf(2,2) = -2*(x(1)^2-1)+x(1)*cos(x(1)*x(2))*(u+1/3*u^3+sin(u));
Gf(4,3) = 1;
Gf(1,4) = 1;
Gf(3,4) = -1;
Gf(4,4) = -0.2;

end