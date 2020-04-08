%plot numerical eigenvalues of space discretization, like in WT2005 paper
%3 branches, etc.

clear all;

%defining globals before calling global_space
global resolution wavenumber mlev dispfig dispom;

%set up parameters to match k value as in TW2005 paper
resolution=1;
wavenumber=360*40/12;

mlev=20;

%plotting
dispfig=20;
dispom=20;

global_space;

thuburn_system_plot_omega();
