% ex 3
% internal parameters
f = 1.378;
kx = 640/1.417;
ky = 480/0.945;

x0 = 320;
y0 = 240;

% Calibration matrix:
K = [f*kx 0 x0; 0 f*ky y0; 0 0 1];

%Kinv = inv(K);

P1x = -0.023;
P1y = -0.261;
P1z = 2.376;

P2x = 0.659;
P2y = -0.071;
P2z = 2.082;

p1x = 52;
p1y = 163;

p2x = 218;
p2y = 216;

syms rho c1 c2 c3;

% P1
p1 = [p1x; p1y; 1];
P1 = [P1x; P1y; P1z];

% P2
p2 = [p2x; p2y; 1];
P2 = [P2x; P2y; P2z];

C = [c1; c2; c3];


C1 = P1 - K\(rho*p1);
C2 = P2 - K\(rho*p2);

% solve(P1 - K\(rho*p1) - C, P2 - K\(rho*p2) - C)

 A = [9047/33920; 63/832];
 b = [341/500; .19];
 rho = A\b;
 
 subs(C1,rho,2.5535)
 subs(C2,rho,2.5535)