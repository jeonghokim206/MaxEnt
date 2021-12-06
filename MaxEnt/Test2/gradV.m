 % Calculate the gradient of V appears in the paper Chow et. al

function dv = gradV(t,x,T)  

% v0=[0;0;0;0];
v0 = 2*rand(4,1)-1;
v=v0;
alpha = 0.03;
epsilon = 0.5*1e-3;
count = 0;
M = 100;

jk=1;
k=1;
while k<=M
    e = zeros(4,1); e(jk)=1;
    vp = v;
    
    gradF = (Fxvt(x,v+1e-3*e,T-t)-Fxvt(x,v,T-t))/1e-3;
    vp(jk) = v(jk)-alpha * gradF;
    jk = jk+1;
    if jk == 5
        jk=1;
    end
    if norm(vp-v)>epsilon
        count = 0;
    end
    if norm(vp-v)<epsilon
        count=count+1;
    end
    v=vp;
    if count==4
        break
    end
    k=k+1;
end

dv = v;

end