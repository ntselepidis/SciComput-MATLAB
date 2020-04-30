clear; clc;
syms xi eta J y41 y12 x14 x21

N1=(1/4)*(1-xi)*(1-eta);
N2=(1/4)*(1+xi)*(1-eta);
N3=(1/4)*(1+xi)*(1+eta);
N4=(1/4)*(1-xi)*(1+eta);

N1x=(1/(8*J))*(y41*(eta-1)+y12*(xi-1));
N1y=(1/(8*J))*(x14*(eta-1)+x21*(xi-1));

N2x=(1/(8*J))*(y41*(-eta+1)+y12*(-xi-1));
N2y=(1/(8*J))*(x14*(-eta+1)+x21*(-xi-1));

N3x=(1/(8*J))*(y41*(eta+1)+y12*(xi+1));
N3y=(1/(8*J))*(x14*(eta+1)+x21*(xi+1));

N4x=(1/(8*J))*(y41*(-eta-1)+y12*(-xi+1));
N4y=(1/(8*J))*(x14*(-eta-1)+x21*(-xi+1));

% Diffusion Terms

W11=J*int(int(N1x*N1x+N1y*N1y,xi,-1,1),eta,-1,1);
W22=J*int(int(N2x*N2x+N2y*N2y,xi,-1,1),eta,-1,1);
W33=J*int(int(N3x*N3x+N3y*N3y,xi,-1,1),eta,-1,1);
W44=J*int(int(N4x*N4x+N4y*N4y,xi,-1,1),eta,-1,1);

W12=J*int(int(N1x*N2x+N1y*N2y,xi,-1,1),eta,-1,1);
W13=J*int(int(N1x*N3x+N1y*N3y,xi,-1,1),eta,-1,1);
W14=J*int(int(N1x*N4x+N1y*N4y,xi,-1,1),eta,-1,1);

W23=J*int(int(N2x*N3x+N2y*N3y,xi,-1,1),eta,-1,1);
W24=J*int(int(N2x*N4x+N2y*N4y,xi,-1,1),eta,-1,1);

W34=J*int(int(N3x*N4x+N3y*N4y,xi,-1,1),eta,-1,1);

% Reaction Terms

M11=J*int(int(N1*N1,xi,-1,1),eta,-1,1);
M22=J*int(int(N2*N2,xi,-1,1),eta,-1,1);
M33=J*int(int(N3*N3,xi,-1,1),eta,-1,1);
M44=J*int(int(N4*N4,xi,-1,1),eta,-1,1);

M12=J*int(int(N1*N2,xi,-1,1),eta,-1,1);
M13=J*int(int(N1*N3,xi,-1,1),eta,-1,1);
M14=J*int(int(N1*N4,xi,-1,1),eta,-1,1);

M23=J*int(int(N2*N3,xi,-1,1),eta,-1,1);
M24=J*int(int(N2*N4,xi,-1,1),eta,-1,1);

M34=J*int(int(N3*N4,xi,-1,1),eta,-1,1);

% Source Terms

f1=J*int(int(N1,xi,-1,1),eta,-1,1);
f2=J*int(int(N2,xi,-1,1),eta,-1,1);
f3=J*int(int(N3,xi,-1,1),eta,-1,1);
f4=J*int(int(N4,xi,-1,1),eta,-1,1);
