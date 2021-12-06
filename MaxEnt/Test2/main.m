%% The main file, computing the controlled/uncontrolled trajectories.


dt = 0.01; % Time step
T = 2.5;   % Maximal time for single iteration
time = 0:dt:T;
iter = length(time);


x_total = zeros(4,8*(iter-1)+1); % prelocating the position vector for the 8 copies of the iteration
x0 = 0.1*[0.5;2.5;0;0.2]; % initial data
x_total(:,1) = x0; % initial data allocation

for K = 1:8
    K
    x0 = x_total(:,(K-1)*(iter-1)+1);
    x = zeros(4,iter);
    u = zeros(1,iter);
    x(:,1) = x0;
    for i = 1:iter -1
        i
%         while true
%             grad = gradV(time(i),x(:,i),T);
%             nancheck = isnan(grad);
%             if sum(nancheck)==0
%                 break
%             end
%         end
        
        % Uncontrolled
        
        U = zeros(1,1000);          
        
        % Controlled
        
%         for j = 1:1000
%            while true
%                 y = 2*rand-1; q=rand(1);
%                  if q<exp(-(dot(grad,vanderpole(x(:,i),y))+norm(x(:,i),1)+abs(y))/1)/exp(2.5/1) || exp(-(dot(grad,vanderpole(x(:,i),y))+norm(x(:,i),1)+abs(y))/1)/exp(2.5/1)<1e-6
%                      U(j)=y;
%                     break
%                  end
%            end
%         end


        force = zeros(4,1000);
        for j = 1:1000
            force(:,j) = vanderpole(x(:,i),U(j));
        end
        u(i) = mean(U);    
        x(:,i+1) = x(:,i) + dt * mean(force,2);

    end
    x_total(:,(K-1)*(iter-1)+2:K*(iter-1)+1) = x(:,2:end);
end