%% colechelon_symbolic_model3.m
%%
%% PURPOSE
%%   Compute quadratic invariants of Model 3 (K=3, M=5)
%%   using the standard approach
%%
%% PAPER REFERENCE
%%   Seshadri & Lakshmivarahan (2026a)
%%
%% MODEL
%%   Three gyrostats sharing modes in a dense pattern (gyrostats 1,2,3
%%   on mode triples that overlap in two modes each). Parameters as above.
%%
%% REQUIREMENTS: MATLAB Symbolic Math Toolbox
%%
%% OUTPUT
%%   A  -- constraint matrix (A8 is its final form after column reduction)

clear, close all

p1 = sym('p1'); q1 = sym('q1'); 
p2 = sym('p2'); q2 = sym('q2');
p3 = sym('p3'); q3 = sym('q3');

r1 = -(p1+q1); r2 = -(p2+q2); 

r3 = -(p3+q3); 

a1 = sym('a1'); b1 = sym('b1'); c1 = sym('c1');
a2 = sym('a2'); b2 = sym('b2'); c2 = sym('c2');
a3 = sym('a3'); b3 = sym('b3'); c3 = sym('c3');


A = [p1 q1 r1 0 0 0 0 0 0 0; ...
    -c1-c3 c1+c3 0 0 0 0 0 r1 r3 0; ...
    0 -a1 a1 0 0 p1 0 0 0 0; ...
    b1 0 -b1 0 0 0 q1 0 0 0; ...
    0 0 p2 q2 r2 0 0 0 0 0; ...
    0 0 -c2 c2 0 0 0 0 0 r2; ...
    0 0 0 -a2 a2 0 0 p2 0 0; ...
    0 0 b2 0 -b2 0 0 0 q2 0; ...
    p3 q3 0 r3 0 0 0 0 0 0; ...
    b3 0 0 -b3 0 0 q3 0 0 0; ...
    0 -a3 0 a3 0 p3 0 0 0 0; ...
    0 0 0 0 0 0 c1+c3 -b1 -b3 0; ...
    0 0 0 0 0 -c1-c3 0 a1 a3 0; ...
    0 0 0 0 0 b1 -a1 0 c2 -b2; ...
    0 0 0 0 0 b3 -a3 -c2 0 a2;... 
    0 0 0 0 0 0 0 b2 -a2 0];

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

alpha8 = - A1(2,8)/A1(2,2);
a8 = simplify(A1(:,8) + alpha8*A1(:,2));

alpha9 = - A1(2,9)/A1(2,2);
a9 = simplify(A1(:,9) + alpha9*A1(:,2));

A2(:,3) = a3; 
A2(:,8) = a8; 
A2(:,9) = a9;

% col. exchange
i1 = 3; i2 = 6;
atemp = A2(:,i1);
A2(:,i1) = A2(:,i2);
A2(:,i2) = atemp; 

%% Third stage
A3 = A2; 

alpha8 = - A2(3,8)/A2(3,3);
a8 = simplify(A2(:,8) + alpha8*A2(:,3));

alpha9 = - A2(3,9)/A2(3,3);
a9 = simplify(A2(:,9) + alpha9*A2(:,3));

A3(:,8) = a8; 
A3(:,9) = a9; 

% col. exchange
i1 = 4; i2 = 7;
atemp = A3(:,i1);
A3(:,i1) = A3(:,i2);
A3(:,i2) = atemp; 

%% Fourth stage
A4 = A3; 

alpha8 = - A3(4,8)/A3(4,4);
a8 = simplify(A3(:,8) + alpha8*A3(:,4));

alpha9 = - A3(4,9)/A3(4,4);
a9 = simplify(A3(:,9) + alpha9*A3(:,4));

A4(:,8) = a8; 
A4(:,9) = a9; 

%% Fifth stage
A5 = A4; 

alpha6 = - A4(5,6)/A4(5,5);
a6 = simplify(A4(:,6) + alpha6*A4(:,5));

alpha7 = - A4(5,7)/A4(5,5);
a7 = simplify(A4(:,7) + alpha7*A4(:,5));

A5(:,6) = a6;
A5(:,7) = a7;

%% Sixth stage
A6 = A5; 

alpha7 = - A5(6,7)/A5(6,6);
a7 = simplify(A5(:,7) + alpha7*A5(:,6));

alpha10 = - A5(6,10)/A5(6,6);
a10 = simplify(A5(:,10) + alpha10*A5(:,6));

A6(:,7) = a7;
A6(:,10) = a10;

% col. exchange
i1 = 7; i2 = 8;
atemp = A6(:,i1);
A6(:,i1) = A6(:,i2);
A6(:,i2) = atemp; 

% col. exchange
i1 = 8; i2 = 10;
atemp = A6(:,i1);
A6(:,i1) = A6(:,i2);
A6(:,i2) = atemp; 

%% Seventh stage
A7 = A6; 

alpha8 = - A6(7,8)/A6(7,7);
a8 = simplify(A6(:,8) + alpha8*A6(:,7));

A7(:,8) = a8;

%% Eighth stage
A8 = A7; 

alpha9 = - A7(8,9)/A7(8,8);
a9 = simplify(simplify(A7(:,9) + simplify(alpha9*A7(:,8))));

A8(:,9) = a9;
