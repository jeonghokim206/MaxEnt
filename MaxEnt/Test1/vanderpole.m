function y = vanderpole(x,u)
    y = zeros(2,1);
    y(1) = x(2);
    y(2) = -2*(x(1)^2-1)*x(2)-x(1)+(2+sin(x(1)*x(2)))*(u+1/3*u^3+sin(u));
end
