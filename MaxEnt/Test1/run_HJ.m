clc, clf, clear, close all

N = 100;
x = linspace(-1,1,N);
y = linspace(-1,1,N);
dx = x(2)-x(1);
dy = y(2)-y(1);
W = zeros(N,N);
for i = 1:N
    for j = 1:N
        W(i,j) = abs(x(i))+abs(y(j));
    end
end

dt = 0.001;
iter = ceil(0.1/dt);
Wp = W;

for k=1:iter-1
    for i = 2:N-1
        for j=2:N-1
            ip = i+1; in = i-1; jp = j+1; jn = j-1;
          
            up = (W(ip,j)-W(i,j))/dx;
            un = (W(i,j)-W(in,j))/dx;
            vp = (W(i,jp)-W(i,j))/dy;
            vn = (W(i,j)-W(i,jn))/dy;
            
%             Wp(i,j) = W(i,j) - dt*HLF(x(i),y(j),up,un,vp,vn);
            Wp(i,j) = W(i,j) - dt*HG(x(i),y(j),up,un,vp,vn);
        end
    end
    for j=2:N-1
%         Wp(1,j) = 2*Wp(2,j)-Wp(3,j);
%         Wp(N,j) = 2*Wp(N-1,j)-Wp(N-2,j);
        Wp(1,j) = 3*Wp(2,j)-3*Wp(3,j)+Wp(4,j);
        Wp(N,j) = 3*Wp(N-1,j)-3*Wp(N-2,j)+Wp(N-3,j);
    end
    for i=2:N-1
%         Wp(i,1) = 2*Wp(i,2)-Wp(i,3);
%         Wp(i,N) = 2*Wp(i,N-1)-Wp(i,N-2);
        Wp(i,1) = 3*Wp(i,2)-3*Wp(i,3)+Wp(i,4);
        Wp(i,N) = 3*Wp(i,N-1)-3*Wp(i,N-2)+Wp(i,N-3);
    end
    Wp(1,1) = Wp(1,2)+Wp(2,1)-Wp(2,2);
    Wp(N,1) = Wp(N-1,1)+Wp(N,2)-Wp(N-1,2);
    Wp(1,N) = Wp(2,N)+Wp(1,N-1)-Wp(2,N-1);
    Wp(N,N) = Wp(N,N-1)+Wp(N-1,N)-Wp(N-1,N-1);
   
    W=Wp;
    k
end
%%
W_grid_free = zeros(N,N);

for i = 1:N
    for j = 1:N
        grad = gradV(0,[x(i);y(j)],0.1);
        W_grid_free(i,j) = Fxvt([x(i);y(j)],grad,0.1);
        i 
        j
    end
end

figure(1)
surf(x,y,W)
zlim([0 2.2])
xlabel('$x$','interpreter','latex','fontsize',18);
ylabel('$y$','interpreter','latex','fontsize',18);
zlabel('$W$','interpreter','latex','fontsize',18);
shading interp
colorbar
caxis([0 2.1])

figure(2)
surf(x,y,W_grid_free)
zlim([0 2.2])
xlabel('$x$','interpreter','latex','fontsize',18);
ylabel('$y$','interpreter','latex','fontsize',18);
zlabel('$W$','interpreter','latex','fontsize',18);
shading interp
colorbar
caxis([0 2.1])


figure(3)
imagesc(x,y,W)
xlabel('$x$','interpreter','latex','fontsize',18);
ylabel('$y$','interpreter','latex','fontsize',18);
colorbar
caxis([0 2.1])


figure(4)
imagesc(x,y,W_grid_free)
xlabel('$x$','interpreter','latex','fontsize',18);
ylabel('$y$','interpreter','latex','fontsize',18);
colorbar
caxis([0 2.1])


figure(5)
imagesc(x,y,W-W_grid_free)
xlabel('$x$','interpreter','latex','fontsize',18);
ylabel('$y$','interpreter','latex','fontsize',18);
colorbar

