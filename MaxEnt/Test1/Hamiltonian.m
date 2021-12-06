% Computing Hamiltonian

function y=Hamiltonian(x,p)

alpha=0.1;

int = 0;
for i = 1:100
     u = -1+0.02*i;
     int = int + exp(-(dot(p,vanderpole(x,u))+(norm(x)+norm(u)))/alpha) * 0.02;
end

y = alpha*log(int);

end