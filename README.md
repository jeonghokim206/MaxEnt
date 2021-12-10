# MaxEnt

This include the codes for the numerical simulation in the paper "Hamilton-Jacobi-Bellman Equations for Maximum Entropy Opitmal Control" in arXiv: https://arxiv.org/abs/2009.13097.
We present codes for three tests:

## Test 1
The first test is to validate the grid-free scheme for Hamilton-Jacobi equations based on the paper "Algorithm for Overcoming the Curse of Dimensionality for State-dependent Hamilton-Jacobi equations" by Y. T. Chow et al. (J. Comp. Phys. 387 pp. 376-409, 2019). We compare the grid-free scheme with the standard finite difference method using Godunov flux. 

In the folder "Test 1", the main file is "run_HJ.m", in which both of the finite difference scheme and grid-free schemes are implemented. Although the numerical results in the arXiv paper is based on the Godunov scheme, we also provide the Lax-Friedrich flux. The files "HG.m" and "HLF.m" are the flux functions for each scheme respectively. On the other hand, the functions "Fxvt.m" and "gradV.m" are the implementation of the grid-free schemes in the paper by Y. T. Chow et al. Finally, the other functions "vanderpole.m" and "Hamiltonian.m" includes the details of the dynamics system and corresponding Hamiltonian.

To conduct the numerical test, just simply run the main file "run_HJ.m". The result returns two matrices "W" and "W_grid_free" which show the numerical solution obtained by the finite difference method and grid-free scheme. Five figures shows the detailed shapes of the numerical solutions, and their difference. On may change the spatial or temporal grid sizes by changing the number of grids "N" and the time step "dt".

---

## Test 2
The second test is to test the maximum entropy optimal control to the four-dimensional extended van der Pol oscillator system, from "Adaptive Control of a Class of Non-Afﬁne Systems using Neural Networks" by B.-J. Yang and A. J. Calise (IEEE Trans. Neural Netw. 18 pp. 1149-1159, 2007). We implement the maximum entropy optimal control by solving the HJB equation, and compare it with the uncontrolled system. The functions "Fxvt.m", "gradV.m", "vanderpole.m" and "Hamiltonian.m" are the same as previous test, while "gradf.m" is used to compute the Jacobian matrix of the dynamics. 

To conduct the numerical test, one only needs to run the main file "main.m". To simulate the uncontrolled system, uncomment the line "U = zeros(1,1000)". Otherwise, to simulate the controlled system, uncomment the lines after the comment "%Controlled". The result returns the trajectory of uncontrolled/controlled system as "x_total". In the folder "Test 2", we present the saved data of the uncontrolled and controlled trajectory as "x_traj_total_uncontrol.mat" and "x_traj_total_control.mat". In order to see the result, one just needs to run the "plot_image.m" file. The results show 1) phase portrait on the (x_1,x_2) spaces, 2) The dynamics of controlled system, and 3) The dynamics of uncontrolled system.

---

## Test 3
The third test is implementation of the maximum entropy optimal control without the model parameters and compare it with the standard sinusoidal exploration. We present two codes for the on-policy and off-policy learning in "linear_ADP_on_policy.m" and "linear_ADP_off_policy.m" respectively. Both of codes are based on the reference "Y. Jiang and Z.-P. Jiang, Robust Adaptive Dynamic Programming, John Wiley & Sons, 2017". Please make sure that the "Control System Toolbox" of the MATLAB is prepared.

To conduct the numerical test, simply run one of the files, depending on the type of test (on-policy or off-policy). One can choose one of the exploration, either MaxEnt exploration or sinusoidal exploration by uncomment each of exploration. The result returns the traj_save that shows the controlled trajectory. Each of main file shows the trajectory of the system before and after the learning the optimal control, and the aggregated total running cost.
