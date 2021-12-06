%% Test

H = zeros(200,200);
x1 = linspace(-10,10,200);
x2 = linspace(-10,10,200);
x3 = 0;
x4 = 0;
p = 20*rand(4,1)-10;
for i = 1:200
    for j=1:200
        H(i,j) = Hamiltonian([x1(i),x2(j),x3,x4],p);
    end
end

surf(H)
colorbar
shading interp

%%

F = zeros(20,20);
v1 = linspace(-20,20,20);
v2 = linspace(-20,20,20);
v3 = 0;
v4 = 0;
x = [-0.8817;0.6523;-0.6203;-0.5505];
for i = 1:20
    for j=1:20
        F(i,j) = Fxvt(x,[v1(i),v2(j),v3,v4],3);
        i
        j
    end
end

imagesc(F)
colorbar
