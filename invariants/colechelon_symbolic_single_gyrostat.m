%% colechelon_symbolic_single_gyrostat.m
%%
%% PURPOSE
%%   Compute the number of quadratic invariants of a single Volterra gyrostat (K=1, M=3)
%%   using the standard approach
%%
%% PAPER REFERENCE
%%   Seshadri & Lakshmivarahan (2026a)
%%
%% METHOD
%%   Builds the 3x7 constraint matrix A whose null space gives the
%%   coefficients of the quadratic invariant C = d_i x_i^2/2 + f_i x_i.
%%   The constraint dC/dt = 0 for all x gives A*u = 0 where
%%   u = (d1,d2,d3, f1,f2,f3, e12). Column echelon reduction then
%%   identifies the free variables and the invariant structure.
%%
%% REQUIREMENTS: MATLAB Symbolic Math Toolbox
%%
%% OUTPUT
%%   A  -- constraint matrix (A4 is its final form after column reduction)

clear, close all

p1 = sym('p1'); q1 = sym('q1'); 

r1 = -(p1+q1);


a1 = sym('a1'); b1 = sym('b1'); c1 = sym('c1');


A = [p1 q1 r1 0 0 0; ...
    -c1 c1 0 0 0 r1; ...
    0 -a1 a1 p1 0 0; ...
    b1 0 -b1 0 q1 0; ...
    0 0 0 0 c1 -b1; ...
    0 0 0 -c1 0 a1; ...
    0 0 0 b1 -a1 0];

%% First stage
A1 = A; 

alpha2 =  - A(1,2)/A(1,1);
a2 = simplify(A(:,2) + alpha2*A(:,1));

alpha3 =  - A(1,3)/A(1,1);
a3 = simplify(A(:,3) + alpha3*A(:,1));

A1(:,2) = a2;
A1(:,3) = a3;

%% Second stage
A2 = A1;

alpha3 =  - A1(2,3)/A1(2,2);
a3 = simplify(A1(:,3) + alpha3*A1(:,2));

alpha6 = - A1(2,6)/A1(2,2);
a6 = simplify(A1(:,6) + alpha6*A1(:,2));

A2(:,3) = a3; 
A2(:,6) = a6; 

% col. exchange
i1 = 3; i2 = 4;
atemp = A2(:,i1);
A2(:,i1) = A2(:,i2);
A2(:,i2) = atemp; 

% col. exchange
i1 = 4; i2 = 6;
atemp = A2(:,i1);
A2(:,i1) = A2(:,i2);
A2(:,i2) = atemp; 

%% Third stage
A3 = A2; 

alpha4 = - A2(3,4)/A2(3,3);
a4 = simplify(A2(:,4) + alpha4*A2(:,3));

A3(:,4) = a4; 

%% Fourth stage
A4 = A3; 

alpha5 = - A3(4,5)/A3(4,4);
a5 = simplify(A3(:,5) + alpha5*A3(:,4));

A4(:,5) = a5; 
