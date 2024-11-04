%{
Main
Authors:
Description: This script applies the finite difference method to
numerically computer derivatives in the Radial Equilibrium Equation.

Inputs: 

Outputs:


%}

% Define the geometry of the axial flow compressor stage. All units in SI

d_hub = 0.9;
d_tip = 1;
alpha_1_h = 32.75*pi/180;
beta_1_h = 55.15*pi/180;
alpha_1_t = 30*pi/180;
beta_1_t = 60*pi/180;
N = 6000*pi/30;
c_x_in = 136;
rho_in = 1.5;
p_in = 100000;
T_in = 25 + 273.15;
loss_coeff = 0.05;

