%% Casimir_gradient_consistency.m
%%
%% PURPOSE
%%   Verify Theorem (Consistency of Casimir gradients) for the sparse
%%   Hamiltonian hierarchy: confirms that the null-space intersection
%%   argument gives a Casimir gradient at level K that is collinear with
%%   the projection of the level K+1 gradient (Appendix 4 of main text).
%%
%% PAPER REFERENCE
%%   Seshadri & Lakshmivarahan (2026a)
%%
%% METHOD
%%   Computes V1 = null(J^(1) extended to R^5) and Vd = null(DeltaJ)
%%   where DeltaJ = J^(K=2) - J^(K=1) is the incremental Poisson matrix.
%%   Solves [V1, -Vd]*[x;y] = 0 for the intersection null vector, then
%%   verifies V1*x equals the known Casimir gradient from Table 4.
%%
%%
%% REQUIREMENTS: MATLAB Symbolic Math Toolbox
%%

clear, close all hidden

% Parameters
p1 = sym('p1'); q1 = sym('q1'); r1 = sym('r1');
p2 = sym('p2'); q2 = 0*sym('q2'); r2 = sym('r2');

a1 = sym('a1'); b1 = sym('b1'); c1 = sym('c1');
a2 = sym('a2'); b2 = sym('b2'); c2 = sym('c2');

% Variables
x1 = sym('x1'); x2 = sym('x2'); x3 = sym('x3'); x4 = sym('x4'); x5 = sym('x5'); x6 = sym('x6'); x7 = sym('x7'); x8 = sym('x8'); x9 = sym('x9');

%% Define each gyrostat's contribution to J

J1 = [0 -c1 p1*x2+b1; ...
    c1 0 q1*x1-a1; ...
    -(p1*x2+b1) -(q1*x1-a1) 0];

J2 = [0 -c2 p2*x4+b2; ...
    c2 0 q2*x3-a2; ...
    -(p2*x4+b2) -(q2*x3-a2) 0];


%% Calculate J and null-spaces


Jb1 = sym(zeros(5,5));
Jb1(1:3,1:3) = Jb1(1:3,1:3) + J1

N1 = null(Jb1)

Jb2 = sym(zeros(5,5));
Jb2(1:3,1:3) = Jb2(1:3,1:3) + J1;
Jb2(3:5,3:5) = Jb2(3:5,3:5) + J2

N2 = null(Jb2)

Jb3 = sym(zeros(5,5));
Jb3(3:5,3:5) = Jb3(3:5,3:5) + J2

N3 = null(Jb3)

dJ = Jb2 - Jb1;
Nd = null(dJ)

%% Intersection of null-spaces

V1 = N1; 
Vd = Nd;

V = cat(2,V1,-Vd)

nV = null(V)

x = nV(1:3,1);
y = nV(4:6,1);

int1 = V1*x
int2 = Vd*y