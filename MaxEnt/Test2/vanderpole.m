% Extended van der Pol oscillator from the Referecne.pdf

function y = vanderpole(x,u)
    y = zeros(4,1);
    y(1) = x(2);
    y(2) = -2*(x(1)^2-1)*x(2)-x(1)+(2+sin(x(1)*x(2)))*(u+1/3*u^3+sin(u));
%     y(2) = -2*(x(1)-1)-x(1)+(u+1/3*u^3+sin(u));
    y(3) = x(4);
    y(4) = -x(3)-0.2*x(4)+x(1);
end
