clear; clc;
syms xi eta J y41 y12 x14 x21

N1=-(1/2)*(xi+eta);
N2=(1/2)*(1+xi);
N4=(1/2)*(1+eta);

N1x=(1/(4*J))*(y41*(-1)+y12*(-1));
N1y=(1/(4*J))*(x14*(-1)+x21*(-1));

N2x=(1/(4*J))*(y41*(1)+y12*(0));
N2y=(1/(4*J))*(x14*(1)+x21*(0));

N4x=(1/(4*J))*(y41*(0)+y12*(1));
N4y=(1/(4*J))*(x14*(0)+x21*(1));

% Diffusion Terms

W11=J*int(int(N1x*N1x+N1y*N1y,xi,-1,-eta),eta,-1,1);
W22=J*int(int(N2x*N2x+N2y*N2y,xi,-1,-eta),eta,-1,1);
W44=J*int(int(N4x*N4x+N4y*N4y,xi,-1,-eta),eta,-1,1);

W12=J*int(int(N1x*N2x+N1y*N2y,xi,-1,-eta),eta,-1,1);
W14=J*int(int(N1x*N4x+N1y*N4y,xi,-1,-eta),eta,-1,1);

W24=J*int(int(N2x*N4x+N2y*N4y,xi,-1,-eta),eta,-1,1);

% Reaction Terms

M11=J*int(int(N1*N1,xi,-1,-eta),eta,-1,1);
M22=J*int(int(N2*N2,xi,-1,-eta),eta,-1,1);
M44=J*int(int(N4*N4,xi,-1,-eta),eta,-1,1);

M12=J*int(int(N1*N2,xi,-1,-eta),eta,-1,1);
M14=J*int(int(N1*N4,xi,-1,-eta),eta,-1,1);

M24=J*int(int(N2*N4,xi,-1,-eta),eta,-1,1);

% Source Terms

f1=J*int(int(N1,xi,-1,-eta),eta,-1,1);
f2=J*int(int(N2,xi,-1,-eta),eta,-1,1);
f4=J*int(int(N4,xi,-1,-eta),eta,-1,1);
