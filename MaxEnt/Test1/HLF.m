function y = HLF(x,y,up,un,vp,vn)
    y = Hamiltonian([x;y],[(up+un)/2;(vp+vn)/2])-0.5*(up-un)-0.5*(vp-vn);
end