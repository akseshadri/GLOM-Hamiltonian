%% Casimir_hierarchy_of_model5.m
%%
%% PURPOSE
%%   Compute Casimir functions and verify Hamiltonian conditions for the
%%   dense GLOM hierarchy motivated by 3D Rayleigh-Benard convection
%%   (Model 5 of the main text).
%%
%% PAPER REFERENCE
%%   Seshadri & Lakshmivarahan (2026a)
%%
%% KEY RESULT
%%   Model 5 is Hamiltonian only under the explicit parameter constraints
%%   of the main text. Adding gyrostats 4 and 5 of this hierarchy eliminates
%%   all Casimirs, illustrating how cross-coupling exhausts the nullspace
%%   of J
%%
%% REQUIREMENTS: MATLAB Symbolic Math Toolbox
%%

clear, close all

p1 = sym('p1'); q1 = 0*sym('q1'); % r1 = sym('r1');
p2 = sym('p2'); q2 = 0*sym('q2'); % r2 = sym('r2');
p3 = 0*sym('p3'); q3 = sym('q3'); % r3 = sym('r3');
p4 = sym('p4'); q4 = 0*sym('q4'); % r4 = sym('r4');
p5 = 0*sym('p5'); q5 = sym('q5'); % r5 = sym('r5');

a1 = sym('a1'); b1 = a1; c1 = a1;
a2 = sym('a2'); b2 = -a2; c2 = -a2;
a3 = sym('a3'); b3 = sym('b3'); c3 = sym('c3');
a4 = sym('a4'); b4 = a4; c4 = -a4;
a5 = sym('a5'); b5 = a5; c5 = a5;

x1 = sym('x1'); x2 = sym('x2'); x3 = sym('x3'); x4 = sym('x4'); x5 = sym('x5'); x6 = sym('x6'); x7 = sym('x7'); x8 = sym('x8');

z = sym('z'); z = 0;

%% Define J

J1 = [0 -c1 p1*x2+b1; ...  % 1 2 3 
    c1 0 q1*x1-a1; ...
    -(p1*x2+b1) -(q1*x1-a1) 0];

J2 = [0 -c2 p2*x4+b2; ...  % 1 4 5
    c2 0 q2*x1-a2; ...
    -(p2*x4+b2) -(q2*x1-a2) 0];

J3 = [0 -c3 p3*x7+b3; ...  % 6 7 8
    c3 0 q3*x6-a3; ...
    -(p3*x7+b3) -(q3*x6-a3) 0];

J4 = [0 -c4 p4*x4+b4; ...  % 3 4 7
    c4 0 q4*x3-a4; ...
    -(p4*x4+b4) -(q4*x3-a4) 0]; 

J5 = [0 -c5 p5*x5+b5; ...  % 2 5 7
    c5 0 q5*x2-a5; ...
    -(p5*x5+b5) -(q5*x2-a5) 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Hierarchy

Jacobiijkmlist = sym(zeros(3,1));

for maincount = 1:5

    if maincount == 5
        J = sym(zeros(8,8));
        J(1:3,1:3) = J1;
        J([1 4 5],[1 4 5]) = J([1 4 5],[1 4 5]) + J2;
        J([6 7 8],[6 7 8]) = J([6 7 8],[6 7 8]) + J3;
        J([3 4 7],[3 4 7]) = J([3 4 7],[3 4 7]) + J4;
        J([2 5 7],[2 5 7]) = J([2 5 7],[2 5 7]) + J5;
        M = 8;

    elseif maincount == 4
        J = sym(zeros(8,8));
        J(1:3,1:3) = J1;
        J([1 4 5],[1 4 5]) = J([1 4 5],[1 4 5]) + J2;
        J([6 7 8],[6 7 8]) = J([6 7 8],[6 7 8]) + J3;
        J([3 4 7],[3 4 7]) = J([3 4 7],[3 4 7]) + J4;
        M = 8;

    elseif maincount == 3
        J = sym(zeros(8,8));
        J(1:3,1:3) = J1;
        J([1 4 5],[1 4 5]) = J([1 4 5],[1 4 5]) + J2;
        J([6 7 8],[6 7 8]) = J([6 7 8],[6 7 8]) + J3;
        M = 8;

    elseif maincount == 2
        J = sym(zeros(5,5));
        J(1:3,1:3) = J1;
        J([1 4 5],[1 4 5]) = J([1 4 5],[1 4 5]) + J2;
        M = 5;

    elseif maincount == 1
        J = sym(zeros(3,3));
        J(1:3,1:3) = J1;
        M = 3;

    end

     %% Levi-Civita matrix

    N = 3;
    [mats{1:N}] = ndgrid(1:N);
    pairsIndex = nchoosek(1:N,2);
    lcMat = sign(prod(cat(N+1,mats{pairsIndex(:,2)})-...
        cat(N+1,mats{pairsIndex(:,1)}),N+1));

    lcMat = sym(lcMat);

    %% Verify that Jacobi identity is satisfied for each member of hierarchy

    Jacobiijkm = 0;

    for i = 1:M
        for j = 1:M
            for k = 1:M

                % alternating tensor
                ijk = [i, j, k];
                uijk = unique(ijk);

                if numel(uijk) <= 2
                    eijk = 0;
                else
                    [sijk,isort] = sort(ijk,'ascend');
                    eijk = lcMat(isort(1),isort(2),isort(3));
                end

                for m = 1:M
                    J_im = J(i,m);

                    J_jk = J(j,k);

                    if m == 1
                        dJ_jk_xm = diff(J_jk,x1);
                    elseif m == 2
                        dJ_jk_xm = diff(J_jk,x2);
                    elseif m == 3
                        dJ_jk_xm = diff(J_jk,x3);
                    elseif m == 4
                        dJ_jk_xm = diff(J_jk,x4);
                    elseif m == 5
                        dJ_jk_xm = diff(J_jk,x5);
                    elseif m == 6
                        dJ_jk_xm = diff(J_jk,x6);
                    elseif m == 7
                        dJ_jk_xm = diff(J_jk,x7);
                    elseif m == 8
                        dJ_jk_xm = diff(J_jk,x8);
                    end
                    Jacobiijkm = Jacobiijkm + eijk*J_im*dJ_jk_xm;

                end
            end
        end
    end

    %% Casimir gradient

    outcount = maincount
    outJid = Jacobiijkm
    outM = M
    NJ = null(J)
    
end


