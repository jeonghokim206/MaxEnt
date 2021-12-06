clc, clf, clear, close all 

% Linear ADP with on-policy learning via Max-ent exploration, compared to
% the sinusoidal exploration. The code is based on the reference "Y. Jiang
% and Z.-P. Jiang, Robust Adaptive Dynamic Programming, John Wiley & Sons,
% 2017"

%% Parameters
xn = 10; % dimensions of state and control
un = 10;

Q = 1e-1*diag(ones(xn,1)); % Q and R matrices
R = 1*diag(ones(un,1));
alpha = 0.1; % relaxation parameter              

rng(200)  % Initialize a fixed model A and B
A=rss(xn,xn,un).A;
A = A + (-max(real(eig(A)))-0.01)*eye(xn);
B=rss(xn,xn,un).B;
A = 1 * A;
B = 0.1 * B;
rng('shuffle')

% Initialize the feedback gain matrix
Kinit  = zeros(un,xn);  %Only if A is Hurwitz, K can be set as zero.
K = Kinit;
i = 0;
MaxIteration = 40;  %Max iteration times
T  = 0.01;          %Length of each integration interval
lambda = 1e-10;

rng(100)
x0 = 40*rand(xn,1)-20; %Initial condition
rng('shuffle')


P_old = zeros(xn);
P = eye(xn)*10;    % Initialize the previous cost matrix
it = 0;            % Counter for iterations

X=[x0;kron(x0,zeros(xn,1));kron(x0,zeros(un,1))];
expl_noise_freq = (rand(un,100)-.5)*100;

data_number_record = [];
traj_save = x0;
cost_save = 0;
t_save = 0;

[Kopt,Popt] = lqr((A-lambda/2*eye(xn)),B,Q,R); % Calculate the ideal solution for comparion purpose

%% On-policy learning
tic
while norm(P-P_old)>5e-1 && it< MaxIteration 
    Theta = [];
    Xi = [];
    
    it = it+1;                        % Update and display the # of iters
    P_old = P;                        % Update the previous cost matrix
    QK = Q+K'*R*K;                    % Update the Qk matrix
    
    x0 = traj_save(:,end);
    X=[x0;kron(x0,zeros(xn,1));kron(x0,zeros(un,1))];
    
    while rank(Theta)< xn*(xn+1)/2 +xn*un 
        
        u = normrnd(-K*X(1:xn),alpha); e = u+K*X(1:xn);                  % MaxEnt exploration
%         u = -Kopt * X(1:xn); e = zeros(un,1);                          % Optimal control
%         e = 0.1*sum(sin(expl_noise_freq*T*i),2); u = -K*X(1:xn)+e;     % Sinusoidal exploration
        
        k1 = A * X(1:xn) + B * u;
        k2 = A * (X(1:xn) + T * k1 / 2) + B * u;
        k3 = A * (X(1:xn) + T * k2 / 2) + B * u;
        k4 = A * (X(1:xn) + T * k3) + B * u;
        xnext = X(1:xn) + T/6 * (k1 + 2*k2 + 2*k3 + k4);
    
        l1 = exp(-lambda*T*i)*kron(X(1:xn),X(1:xn));
        l2 = exp(-lambda*T*(i+0.5))*kron(X(1:xn) + T * k1/2 , X(1:xn) + T * k1/2);
        l3 = exp(-lambda*T*(i+0.5))*kron(X(1:xn) + T * k2/2 , X(1:xn) + T * k2/2);
        l4 = exp(-lambda*T*(i+1))*kron(X(1:xn) + T * k3   , X(1:xn) + T * k3  );
        dxx = X(xn+1:xn+xn^2) + T/6 * (l1 + 2*l2 + 2*l3 + l4);

        m1 = exp(-lambda*T*i)*kron(X(1:xn)            , e);
        m2 = exp(-lambda*T*(i+0.5))*kron(X(1:xn) + T * k1/2 , e);
        m3 = exp(-lambda*T*(i+0.5))*kron(X(1:xn) + T * k2/2 , e);
        m4 = exp(-lambda*T*(i+1))*kron(X(1:xn) + T * k3   , e);
        dxe = X(xn+xn^2+1:end) + T/6 * (m1 + 2*m2 + 2*m3 + m4);
        
        Theta = [Theta ; [exp(-lambda*T*(i+1))*kron(xnext,xnext)'-exp(-lambda*T*i)*kron(X(1:xn),X(1:xn))',-2*(dxe-X(xn+xn^2+1:end))']];
        Xi    = [Xi    ; -(dxx-X(xn+1:xn+xn^2))'*QK(:)];
        X=[xnext;dxx;dxe];
        traj_save = [traj_save,X(1:xn)];
        cost_save = [cost_save,cost_save(end) + exp(-lambda*T*(sum(data_number_record)+i)) * T * (xnext'*Q*xnext + u'*R*u)]; 
        t_save = [t_save, t_save(end)+T];
        i = i +1;
    end
    data_number_record = [data_number_record, i];
    i=0;

    
    pp = pinv(Theta)*Xi;             % Solve the equations in the LS sense
    P = reshape(pp(1:xn*xn), [xn, xn]);  % Reconstruct the symmetric matrix
    P = (P + P')/2;
    
    BPv = pp(end-(xn*un-1):end);
    K = reshape(BPv,un,xn);% Get the improved gain matrix
    
    disp(['K_', num2str(it), '=']);
    disp(K);
    disp(norm(K-Kopt));

end
toc

%% Use the control

t_end = 500;
T_use = 0.01;
t = t_save(end):T_use:t_end;
x = traj_save(:,end);
cost = cost_save(end);
for i = 1:length(t)-1
    
%     u = normrnd(-K*x,alpha); % MaxEnt control
    u = -K*x;
%     u = -Kopt*x;  % Optimal control
%     u = zeros(un,1); % Uncontrolled
    
    k1 = A * x + B * u;
    k2 = A * (x + T_use * k1 / 2) + B * u;
    k3 = A * (x + T_use * k2 / 2) + B * u;
    k4 = A * (x + T_use * k3) + B * u;
    xnext = x + T_use/6 * (k1 + 2*k2 + 2*k3 + k4);    
    traj_save=[traj_save,xnext];
    cost_save = [cost_save,cost_save(end) + exp(-lambda*T*(sum(data_number_record)+i))*T_use * (xnext'*Q*xnext + u'*R*u)]; 
    t_save = [t_save,t_save(end)+T_use];
    x = xnext;
end
cost = cost_save(end);

for i = 1 : length(t_save)
    if max(max(abs(traj_save(:,i:end))))<1
        settling_time = t_save(i);
        break
    end
end

updated_time = [];
for i =1 : length(data_number_record)
    updated_time = [updated_time,sum(data_number_record(1:i))*T];
end

figure(1)
plot(t_save,traj_save,'linewidth',1.5)
xlim([0 20])
ylim([-20 15])
xlabel('$t$','interpreter','latex','fontsize',18)
ylabel('$x$','interpreter','latex','fontsize',18)
hold on 
for i = length(data_number_record)
    line([updated_time(i),updated_time(i)],[-20 15],'color','black','linestyle','--','linewidth',0.25);
end
grid on
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
xlabel('$t$','interpreter','latex','fontsize',24)
ylabel('$x$','interpreter','latex','fontsize',24)
grid on
daspect([1 1.75 1])

figure(2)
plot(t_save,cost_save,'linewidth',1.5)
xlabel('$t$','interpreter','latex','fontsize',18)
ylabel('cost','fontsize',18)
% hold on 
% for i = [1,length(data_number_record)]
%     line([updated_time(i),updated_time(i)],[0,120],'color','black','linestyle','--','linewidth',0.25);
% end
xlim([0 20])
grid on

%% Display results
disp('Approximate Cost Matrix')
P
disp('Optimal Cost Matrix')
Popt
disp('Approximate Gain Matrix')
K
disp('Optimal Gain Matrix')
Kopt
disp('diff')
norm(K-Kopt,'fro')
cost_save(end)
sum(data_number_record)/length(data_number_record)
