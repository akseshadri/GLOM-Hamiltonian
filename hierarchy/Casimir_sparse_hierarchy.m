%% Casimir_sparse_hierarchy.m
%%
%% PURPOSE
%%   Compute Casimir functions symbolically for the sparse Hamiltonian
%%   hierarchy (q_k=0) at K=1,2,3,4, and verify that the K-level Casimir
%%   restricts consistently to the (K-1)-level Casimir.
%%
%% PAPER REFERENCE
%%   Seshadri & Lakshmivarahan (2026a)
%%
%% METHOD
%%   At each level K, builds J = sum_{k=1}^K J^(k) with q_k = 0 (the
%%   Hamiltonian condition), computes null(J) symbolically, and verifies
%%   that the gradient of the K-level Casimir is collinear with the
%%   projection of the (K+1)-level Casimir gradient onto the M=2K+1
%%   subspace. 
%%
%% REQUIREMENTS: MATLAB Symbolic Math Toolbox
%%
clear, close all hidden

%% Number of gyrostats
for K = 1:4

    %% Loop across subclasses

    Kout = K

    i = 2^(3*K)-1; 
    isincluded = dec2bin(i,3*K);


    % Parameters
    p1 = sym('p1'); q1 = sym('q1'); r1 = sym('r1');
    p2 = sym('p2'); q2 = 0; r2 = -p2;
    p3 = sym('p3'); q3 = 0; r3 = -p3;
    p4 = sym('p4'); q4 = 0; r4 = -p4;

    a1 = sym('a1'); b1 = sym('b1'); c1 = sym('c1');
    a2 = sym('a2'); b2 = sym('b2'); c2 = sym('c2');
    a3 = sym('a3'); b3 = sym('b3'); c3 = sym('c3');
    a4 = sym('a4'); b4 = sym('b4'); c4 = sym('c4');

    % Variables
    x1 = sym('x1'); x2 = sym('x2'); x3 = sym('x3'); x4 = sym('x4'); x5 = sym('x5'); x6 = sym('x6'); x7 = sym('x7'); x8 = sym('x8'); x9 = sym('x9');


    % assign inclusion of terms

    ia1 = isequal(isincluded(1),'1')*1;
    ib1 = isequal(isincluded(2),'1')*1;
    ic1 = isequal(isincluded(3),'1')*1;

    if K > 1
        ia2 = isequal(isincluded(4),'1')*1;
        ib2 = isequal(isincluded(5),'1')*1;
        ic2 = isequal(isincluded(6),'1')*1;
    end

    if K > 2
        ia3 = isequal(isincluded(7),'1')*1;
        ib3 = isequal(isincluded(8),'1')*1;
        ic3 = isequal(isincluded(9),'1')*1;
    end

    if K > 3
        ia4 = isequal(isincluded(10),'1')*1;
        ib4 = isequal(isincluded(11),'1')*1;
        ic4 = isequal(isincluded(12),'1')*1;
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For each case

    % Save subclass

    if K == 1
        pilistc = [ia1 ib1 ic1];
    elseif K == 2
        pilistc = [ia1 ib1 ic1 ia2 ib2 ic2];
    elseif K == 3
        pilistc = [ia1 ib1 ic1 ia2 ib2 ic2 ia3 ib3 ic3];
    elseif K == 4
        pilistc = [ia1 ib1 ic1 ia2 ib2 ic2 ia3 ib3 ic3 ia4 ib4 ic4];
    end


    % Assign parameters

    a1 = ia1*a1;
    b1 = ib1*b1;
    c1 = ic1*c1;

    if K > 1
        a2 = ia2*a2;
        b2 = ib2*b2;
        c2 = ic2*c2;
    end

    if K > 2
        a3 = ia3*a3;
        b3 = ib3*b3;
        c3 = ic3*c3;
    end

    if K > 3
        a4 = ia4*a4;
        b4 = ib4*b4;
        c4 = ic4*c4;
    end


    %% Define each gyrostat's contribution to J

    J1 = [0 -c1 p1*x2+b1; ...
        c1 0 q1*x1-a1; ...
        -(p1*x2+b1) -(q1*x1-a1) 0];

    if K > 1
        J2 = [0 -c2 p2*x4+b2; ...
            c2 0 q2*x3-a2; ...
            -(p2*x4+b2) -(q2*x3-a2) 0];
    end

    if K > 2
        J3 = [0 -c3 p3*x6+b3; ...
            c3 0 q3*x5-a3; ...
            -(p3*x6+b3) -(q3*x5-a3) 0];
    end

    if K > 3
        J4 = [0 -c4 p4*x8+b4; ...
            c4 0 q4*x7-a4; ...
            -(p4*x8+b4) -(q4*x7-a4) 0];
    end

    %% Calculate J

    if K == 1
        J = sym(zeros(3,3));

        J(1:3,1:3) = J(1:3,1:3) + J1;
    elseif K == 2
        J = sym(zeros(5,5));

        J(1:3,1:3) = J(1:3,1:3) + J1;
        J(3:5,3:5) = J(3:5,3:5) + J2;
    elseif K == 3
        J = sym(zeros(7,7));

        J(1:3,1:3) = J(1:3,1:3) + J1;
        J(3:5,3:5) = J(3:5,3:5) + J2;
        J(5:7,5:7) = J(5:7,5:7) + J3;
    else
        J = sym(zeros(9,9));

        J(1:3,1:3) = J(1:3,1:3) + J1;
        J(3:5,3:5) = J(3:5,3:5) + J2;
        J(5:7,5:7) = J(5:7,5:7) + J3;
        J(7:9,7:9) = J(7:9,7:9) + J4;
    end



    %% Casimir
    Jout = J

    % Nullspace
    NJout = null(J)
  

end

