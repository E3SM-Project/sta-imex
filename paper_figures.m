
%This script produces figures for GMDD submission, 
%in the same order as in the paper. 
%For dispersion/dissipation plots, download 
%Cao, Y.: Munkres Assignment Algorithm, https://www.mathworks.com/matlabcentral/fileexchange/20328-munkres-assignment-algorithm
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 1, spectrum of system (1-5), 3 branches of waves, 
%Rossby, gravity, acoustics.
%Do not comment out this run if running for the 1st time,
%it is needed for dispersion runs.
plot_omega



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 2a, M1 method, acoustic system with (kx,kz)*dt axes
clear all; global which; which=10030;
%high res plot as in paper
%nnx=400; nnz=400;
%low res plot
nnx=40; nnz=40;
lock_plotting_courant_axes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 2b, M1 method, acoustic system with (T_x,dt) axes
clear all; global which; which=10030;
%high res plot as in paper
%nnx=100; nnz=100;
%low res plot
nnx=20; nnz=20;
lock_plotting_new

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 2c, M1 method, system of normal modes with (T_x,dt) axes
clear all; global which; which=10030;
%high res plot as in paper
%nnx=100; nnz=100;
%low res plot
nnx=20; nnz=20;
plot_stability



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 3a, ARK2(2,3,2) method, system of normal modes with (T_x,dt) axes
clear all; global which; which=10232;
%high res plot as in paper
%nnx=100; nnz=100;
%low res plot
nnx=20; nnz=20;
plot_stability

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 3b, IMEX-SSP2(2,3,2) method, system of normal modes with (T_x,dt) axes
clear all; global which; which=10532;
%high res plot as in paper
%nnx=100; nnz=100;
%low res plot
nnx=20; nnz=20;
plot_stability



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 4a, M2a method, system of normal modes with (T_x,dt) axes
clear all; global which; which=10100;
%high res plot as in paper
%nnx=100; nnz=100;
%low res plot
nnx=20; nnz=20;
plot_stability

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 4b, M2b method, system of normal modes with (T_x,dt) axes
clear all; global which; which=10101;
%high res plot as in paper
%nnx=100; nnz=100;
%low res plot
nnx=20; nnz=20;
plot_stability

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 4c, M2c method, system of normal modes with (T_x,dt) axes
clear all; global which; which=10102;
%high res plot as in paper
%nnx=100; nnz=100;
%low res plot
nnx=20; nnz=20;
plot_stability

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 4d, M2a method, dispersion/dissipation
clear all; global which; which=10100;
dispersion_dissipation4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 4e, M2b method, dispersion/dissipation
clear all; global which; which=10101;
dispersion_dissipation4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 4f, M2c method, dispersion/dissipation
clear all; global which; which=10102;
dispersion_dissipation4



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 5a, M2be method, system of normal modes with (T_x,dt) axes
clear all; global which; which=10103;
%high res plot as in paper
%nnx=100; nnz=100;
%low res plot
nnx=20; nnz=20;
plot_stability

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 5b, M2cn method, system of normal modes with (T_x,dt) axes
clear all; global which; which=10104;
%high res plot as in paper
%nnx=100; nnz=100;
%low res plot
nnx=20; nnz=20;
plot_stability

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 5c, M2cno method, system of normal modes with (T_x,dt) axes
clear all; global which; which=10105;
%high res plot as in paper
%nnx=100; nnz=100;
%low res plot
nnx=20; nnz=20;
plot_stability

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 5d, M2be method, dispersion/dissipation
clear all; global which; which=10103;
dispersion_dissipation4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 5e, M2cn method, dispersion/dissipation
clear all; global which; which=10104;
dispersion_dissipation4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 5f, M2cno method, dispersion/dissipation
clear all; global which; which=10105;
dispersion_dissipation4



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 6a, M1 method, convergence wrt dz
clear all; global which; which=10030;
%high res plot as in paper
%nndt=100; nndz=50;
%low res plot
nndt=20; nndz=20;
nlev_convergence

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 6b, ARK2(2,3,2) method, convergence wrt dz
clear all; global which; which=10232;
%high res plot as in paper
%nndt=100; nndz=50;
%low res plot
nndt=20; nndz=20;
nlev_convergence

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figure 6c, M2b method, convergence wrt dz
clear all; global which; which=10101;
%high res plot as in paper
%nndt=100; nndz=50;
%low res plot
nndt=20; nndz=20;
nlev_convergence
