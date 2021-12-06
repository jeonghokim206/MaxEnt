function y = HG(x,y,up,un,vp,vn)
%     H1 = Hamiltonian([x;y],[up;vp]);
%     H2 = Hamiltonian([x;y],[up;vn]);
%     H3 = Hamiltonian([x;y],[un;vp]);
%     H4 = Hamiltonian([x;y],[un;vn]);
%     
%     if un<up && vn<vp
%         y = min(min(H1,H2),min(H3,H4));
%     elseif un<up && vp<vn
%         y = min(max(H1,H2),max(H3,H4));
%     elseif up<un && vn<vp
%         y = max(min(H1,H2),min(H3,H4));
%     else
%         y = max(max(H1,H2),max(H3,H4));
%     end
    H1 = Hamiltonian([x;y],[un;vn]);
    H2 = Hamiltonian([x;y],[un;(vn+vp)/2]);
    H3 = Hamiltonian([x;y],[un;vp]);
    H4 = Hamiltonian([x;y],[(un+up)/2;vn]);
    H5 = Hamiltonian([x;y],[(un+up)/2;(vn+vp)/2]);
    H6 = Hamiltonian([x;y],[(un+up)/2;vp]);
    H7 = Hamiltonian([x;y],[up;vn]);
    H8 = Hamiltonian([x;y],[up;(vn+vp)/2]);
    H9 = Hamiltonian([x;y],[up;vp]);
    
    if un<up && vn<vp
        y = min([min([H1,H2,H3]),min([H4,H5,H6]),min([H7,H8,H9])]);
    elseif un<up && vp<vn
        y = min([max([H1,H2,H3]),max([H4,H5,H6]),max([H7,H8,H9])]);
    elseif up<un && vn<vp
        y = max([min([H1,H2,H3]),min([H4,H5,H6]),min([H7,H8,H9])]);
    else
        y = max([max([H1,H2,H3]),max([H4,H5,H6]),max([H7,H8,H9])]);
    end 

end