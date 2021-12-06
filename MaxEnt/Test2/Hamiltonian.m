function y=Hamiltonian(x,p)

alpha=1;

U = zeros(100,1);
for i = 1:100
     u = -1+0.02*i;
     U(i) = -(dot(p,vanderpole(x,u))+(norm(x,1)+norm(u,1)));
end
K = max(U);
U = U-K;
y = alpha*(K+log(sum(exp(U))*0.02));
% 
% int = 0;
% for i = 1:100
%      u = -1+0.02*i;
%      int = int + exp(-(dot(p,vanderpole(x,u))+(norm(x,1)+norm(u,1)))/alpha) * 0.02;
% end
% 
% y = alpha*log(int);

end