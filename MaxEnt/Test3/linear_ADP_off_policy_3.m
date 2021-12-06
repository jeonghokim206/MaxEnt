clc, clf, clear, close all % 2020-09-14

%% Parameters configuration
xn = 20;
un = 20;

% Set the weighting matrices for the cost function
Q = 0.1*diag(ones(xn,1));
R = diag(ones(un,1));
alpha = 1;

rng(150) % Initialize a fixed model A and B
A=rss(xn,xn,un).A;
A = A + (-max(real(eig(A)))-0.01)*eye(xn);
B=rss(xn,xn,un).B;
A = A;
B = 0.1*B;
rng('shuffle')

% Initialize the feedback gain matrix
Kinit  = zeros(un,xn);  %Only if A is Hurwitz, K can be set as zero.
K = Kinit;
MaxIteration = 40;  %Max iteration times
T  = 0.01;          %Length of each integration interval
i = 0;

lambda = 1e-10;

rng(200)
x0 = 40*rand(xn,1)-20; %Initial condition
rng('shuffle')


% Matrices to collect online data and perform learning
Dxx=[];
Ixx=[];
Ixu=[];

traj_save = x0;
cost_save = 0;
t_save = 0;

% Initial condition of the augmented system
X=[x0;kron(x0,zeros(xn,1));kron(x0,zeros(un,1))];
expl_noise_freq = (rand(un,100)-.5)*100;
%% Run the simulation and obtain the data matrices \delta_{xx},
%I_{xx}, and I_{xu}

tic
while rank([Ixx, Ixu])<xn*(xn+1)/2 +xn*un
    
%     Simulation the system and at the same time collect online info.
    u = normrnd(-K*X(1:xn),alpha);
%     u = 0.1*sum(sin(expl_noise_freq*T*i),2);
    
    k1 = A * X(1:xn) + B * u;
    k2 = A * (X(1:xn) + T * k1 / 2) + B * u;
    k3 = A * (X(1:xn) + T * k2 / 2) + B * u;
    k4 = A * (X(1:xn) + T * k3) + B * u;
    xnext = X(1:xn) + T/6 * (k1 + 2*k2 + 2*k3 + k4);
    
    l1 = exp(-lambda*T*i)*kron(X(1:xn),X(1:xn));
    l2 = exp(-lambda*T*(i+0.5))*kron(X(1:xn)+ T*k1/2,X(1:xn)+ T*k1/2);
    l3 = exp(-lambda*T*(i+0.5))*kron(X(1:xn)+ T*k2/2,X(1:xn)+ T*k2/2);
    l4 = exp(-lambda*T*(i+1))*kron(X(1:xn)+ T*k3,X(1:xn)+ T*k3);
    dxx = X(xn+1:xn+xn^2) + T/6 * (l1 + 2*l2 + 2*l3 + l4);
    
    m1 = exp(-lambda*T*i)*kron(X(1:xn),u);
    m2 = exp(-lambda*T*(i+0.5))*kron(X(1:xn)+ T*k1/2,u);
    m3 = exp(-lambda*T*(i+0.5))*kron(X(1:xn)+ T*k2/2,u);
    m4 = exp(-lambda*T*(i+1))*kron(X(1:xn)+ T*k3,u);
    dxu = X(xn+xn^2+1:end) + T/6 * (m1 + 2*m2 + 2*m3 + m4);

    
    %Append new data to the data matrices
    Dxx=[Dxx;exp(-lambda*T*(i+1))*kron(xnext,xnext)'-exp(-lambda*T*i)*kron(X(1:xn),X(1:xn))'];
    Ixx=[Ixx;dxx'-X(xn+1:xn+xn^2)'];
    Ixu=[Ixu;dxu'-X(xn+xn^2+1:end)'];
    
    X=[xnext;dxx;dxu];
    traj_save = [traj_save,X(1:xn)];
    cost_save = [cost_save,cost_save(end) + exp(-lambda*T*i) * T * (xnext'*Q*xnext + u'*R*u)]; 
    t_save = [t_save, t_save(end)+T];
    i = i+1;
    if mod(i,100)==0
        i
    end
end
N = i;

P_old = zeros(xn);
P = eye(xn)*10;    % Initialize the previous cost matrix
it = 0;            % Counter for iterations
p_save = [];       % Track the cost matrices in all the iterations
k_save = [];       % Track the feedback gain matrix in each iterations

% The system matrices are copied here only for the purpose of analyzing the
% results


[Kopt,Popt] = lqr((A-lambda/2*eye(xn)),B,Q,R); % Calculate the ideal solution for comparion purpose
k_save  = norm(K-Kopt);  % keep track of the differences between the actual K
% and the idea valu

%% Off-policy learning using the collected online data
while norm(P-P_old)>1e-3 & it<MaxIteration
    it = it+1;                        % Update and display the # of iters
    P_old = P;                        % Update the previous cost matrix
    QK = Q+K'*R*K;                    % Update the Qk matrix
    
    Theta = [Dxx,-Ixx*kron(eye(xn),K')-Ixu];  % Left-hand side of
    % the key equation
    Xi = -Ixx*QK(:);                 % Right-hand side of the key equation
    pp = pinv(Theta)*Xi;             % Solve the equations in the LS sense
    P = reshape(pp(1:xn*xn), [xn, xn]);  % Reconstruct the symmetric matrix
    P = (P + P')/2;
    
    BPv = pp(end-(xn*un-1):end);
    K = inv(R)*reshape(BPv,un,xn)/2;% Get the improved gain matrix
    p_save = [p_save,norm(P-Popt)];   % Keep track of the cost matrix
    k_save = [k_save,norm(K-Kopt)];    % Keep track of the control gains
    
    disp(['K_', num2str(it), '=']);
    disp(K);

end
toc
%% Use the control

t_end = 500;
T_test = 0.01;
t = t_save(end):T_test:t_end;
x = traj_save(:,end);

for i = 1:length(t)-1
%     u = -K*x;
    u = -Kopt*x;
%     u = zeros(un,1);
    k1 = A * x + B * u;
    k2 = A * (x + T * k1 / 2) + B * u;
    k3 = A * (x + T * k2 / 2) + B * u;
    k4 = A * (x + T * k3) + B * u;
    xnext = x + T/6 * (k1 + 2*k2 + 2*k3 + k4);        
    traj_save = [traj_save,xnext];
    cost_save = [cost_save,cost_save(end) + exp(-lambda*T*(N+i)) * T_test * (xnext'*Q*xnext + u'*R*u)]; 
    t_save = [t_save,t_save(end)+T_test];
    x=xnext;
end


for i = 1 : length(t_save)
    if max(max(traj_save(:,i:end)))<1
        settling_time = t_save(i);
        break
    end
end

figure(1)
plot(t_save,traj_save,'linewidth',1.5)
ylim([-30 25])
xlim([0 80])
line([N*T,N*T],[-30 25],'color','black','linestyle','--','linewidth',0.25);
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
xlabel('$t$','interpreter','latex','fontsize',24)
ylabel('$x$','interpreter','latex','fontsize',24)
grid on
daspect([1 0.7 1])

figure(2)
plot(t_save,cost_save)
xlabel('$t$','interpreter','latex','fontsize',18)
ylabel('cost','fontsize',18)
xlim([0 10])
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