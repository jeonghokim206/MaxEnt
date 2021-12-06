% Computing F(x,v,t) appears in the paper of Chow et al.

function y = Fxvt(x,v,t)
    dt = 0.05;
    T = 0:dt:t;
    iter = length(T);
    gamma = zeros(4,iter);
    p = zeros(4,iter);
    gamma(:,end) = x; p(:,end)=v;
    for i = 1:iter-1
        gradx=zeros(4,1);
        gradp=zeros(4,1);
        for j =1:4
            ej = zeros(4,1); ej(j) = 1;
            gradx(j) = (Hamiltonian(gamma(:,end-i+1)+1e-3*ej,p(:,end-i+1))-Hamiltonian(gamma(:,end-i+1),p(:,end-i+1)))/1e-3;
            gradp(j) = (Hamiltonian(gamma(:,end-i+1),p(:,end-i+1)+1e-3*ej)-Hamiltonian(gamma(:,end-i+1),p(:,end-i+1)))/1e-3;
        end
        gamma(:,end-i) = gamma(:,end-i+1)-dt*gradp;
        p(:,end-i) = p(:,end-i+1)+dt*gradx;
%         int1 = 0;
%         int2 = zeros(4,1);
%         int3 = zeros(4,1);
%         for j = 1:100
%             u = -1+0.02*j;
%             int1 = int1 + exp(-(dot(p(:,end-i+1),vanderpole(gamma(:,end-i+1),u))+(norm(gamma(:,end-i+1),1)+norm(u)))/1) * 0.02;
%             int2 = int2 + vanderpole(gamma(:,end-i+1),u)*exp(-(dot(p(:,end-i+1),vanderpole(gamma(:,end-i+1),u))+(norm(gamma(:,end-i+1),1)+norm(u)))/1) * 0.02;
%             int3 = int3 + (gradf(gamma(:,end-i+1),u)'*p(:,end-i+1)+0.2*gamma(:,end-i+1))*exp(-(dot(p(:,end-i+1),vanderpole(gamma(:,end-i+1),u))+(norm(gamma(:,end-i+1),1)+norm(u)))/1) * 0.02;
%         end
%         gradx=-int3/int1;
%         gradp=-int2/int1;
%         
%         gamma(:,end-i) = gamma(:,end-i+1)-dt*gradp;
%         p(:,end-i) = p(:,end-i+1)+dt*gradx;
    end
    
    y = norm(gamma(:,1),1);
    for i =1 : iter-1
        y = y + dt* (dot(p(:,i),(gamma(:,i+1)-gamma(:,i))/dt)-Hamiltonian(gamma(:,i),p(:,i)));
    end

end