%% Casimir_hierarchy_of_model4.m
%%
%% PURPOSE
%%   Compute Casimir functions and verify Hamiltonian conditions for the
%%   GLOM hierarchy motivated by 2D Rayleigh-Benard convection
%%   (Model 4 of the main text), at K=1,2,3 gyrostats.
%%
%% PAPER REFERENCE
%%   Seshadri & Lakshmivarahan (2026a)
%%
%% MODEL
%%   Gluhovsky-type 2D RBC model. Gyrostats couple in a hub-spoke pattern
%%   through the large-scale shear mode x1. Parameters correspond to
%%   buoyancy-shear coupling strengths.
%%
%% KEY RESULT
%%   Model 4 is not Hamiltonian in its natural form (Gluhovsky 2006) without additional parameter constraints.
%%
%% REQUIREMENTS: MATLAB Symbolic Math Toolbox
%%

clear, close all

p1 = sym('p1'); q1 = sym('q1'); % r1 = sym('r1');
p2 = 0*sym('p2'); q2 = 0*sym('q2'); % r2 = sym('r2');
p3 = 0*sym('p3'); q3 = 0*sym('q3'); % r3 = sym('r3');

a1 = sym('a1'); b1 = sym('b1'); c1 = sym('c1');
a2 = sym('a2'); b2 = sym('b2'); c2 = b2;
a3 = sym('a3'); b3 = sym('b3'); c3 = b3;

x1 = sym('x1'); x2 = sym('x2'); x3 = sym('x3'); x4 = sym('x4'); x5 = sym('x5'); x6 = sym('x6'); x7 = sym('x7');

z = sym('z'); z = 0;

%% Define J

J1 = [0 -c1 p1*x2+b1; ...
    c1 0 q1*x1-a1; ...
    -(p1*x2+b1) -(q1*x1-a1) 0];

J2 = [0 -c2 p2*x4+b2; ...
    c2 0 q2*x1-a2; ...
    -(p2*x4+b2) -(q2*x1-a2) 0];

J3 = [0 -c3 p3*x6+b3; ...
    c3 0 q3*x1-a3; ...
    -(p3*x6+b3) -(q3*x1-a3) 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Hierarchy

Jacobiijkmlist = sym(zeros(3,1));

for maincount = 1:1

    if maincount == 1
        J = sym(zeros(7,7));
        J(1:3,1:3) = J1;
        J([1 4 5],[1 4 5]) = J([1 4 5],[1 4 5]) + J2;
        J([1 6 7],[1 6 7]) = J([1 6 7],[1 6 7]) + J3;
        M = 7;

    elseif maincount == 2
        J = sym(zeros(5,5));
        J(1:3,1:3) = J1;
        J([1 4 5],[1 4 5]) = J([1 4 5],[1 4 5]) + J2;
        M = 5;
    else
        J = sym(zeros(3,3));
        J(1:3,1:3) = J1;
        M = 3;
    end

    
    %% Casimir gradient
    Mout = M
    Nullout = null(J)

end


