%% Jacobi_dense_hierarchy.m
%%
%% PURPOSE
%%   Verify the Jacobi condition for the canonical (L_A) representation
%%   of the dense (non-sparse) GLOM hierarchy at K=1,...,6.
%%
%% PAPER REFERENCE
%%   Seshadri & Lakshmivarahan (2026a)
%%
%% MODEL
%%   Dense coupling: gyrostats share two modes (e.g. {x1,x2,x3} and
%%   {x2,x3,x4}), unlike the sparse hierarchy where only one mode is shared.
%%
%%
%% REQUIREMENTS: MATLAB Symbolic Math Toolbox
%%
clear, close all

p1 = sym('p1'); q1 = sym('q1'); 
p2 = sym('p2'); q2 = sym('q2'); 
p3 = sym('p3'); q3 = sym('q3');
p4 = sym('p4'); q4 = sym('q4'); 
p5 = sym('p5'); q5 = sym('q5'); 
p6 = sym('p6'); q6 = sym('q6');

a1 = sym('a1'); b1 = sym('b1'); c1 = sym('c1');
a2 = sym('a2'); b2 = sym('b2'); c2 = sym('c2');
a3 = sym('a3'); b3 = sym('b3'); c3 = sym('c3');
a4 = sym('a4'); b4 = sym('b4'); c4 = sym('c4');
a5 = sym('a5'); b5 = sym('b5'); c5 = sym('c5');
a6 = sym('a6'); b6 = sym('b6'); c6 = sym('c6');

%% State variables
x1 = sym('x1'); x2 = sym('x2'); x3 = sym('x3'); x4 = sym('x4'); x5 = sym('x5');x6 = sym('x6'); x7 = sym('x7'); x8 = sym('x8'); x9 = sym('x9');

z = sym('z'); z = 0;

%% Define each gyrostat's contribution to J

J1 = [0 -c1 p1*x2+b1; ...
    c1 0 q1*x1-a1; ...
    -(p1*x2+b1) -(q1*x1-a1) 0];

J2 = [0 -c2 p2*x3+b2; ...
    c2 0 q2*x2-a2; ...
    -(p2*x3+b2) -(q2*x2-a2) 0];

J3 = [0 -c3 p3*x4+b3; ...
    c3 0 q3*x3-a3; ...
    -(p3*x4+b3) -(q3*x3-a3) 0];

J4 = [0 -c4 p4*x5+b4; ...
    c4 0 q4*x4-a4; ...
    -(p4*x5+b4) -(q4*x4-a4) 0];

J5 = [0 -c5 p5*x6+b5; ...
    c5 0 q5*x5-a5; ...
    -(p5*x6+b5) -(q5*x5-a5) 0];

J6 = [0 -c6 p6*x7+b6; ...
    c6 0 q6*x6-a6; ...
    -(p6*x7+b6) -(q6*x6-a6) 0];

%% For each K

Jacobiijkmlist = sym(zeros(6,1));

for K = 1:6

    outcount = K


    if K == 1
        J = sym(zeros(3,3));

        J(1:3,1:3) = J(1:3,1:3) + J1;
    elseif K == 2
        J = sym(zeros(4,4));

        J(1:3,1:3) = J(1:3,1:3) + J1;
        J(2:4,2:4) = J(2:4,2:4) + J2;
    elseif K == 3
        J = sym(zeros(5,5));

        J(1:3,1:3) = J(1:3,1:3) + J1;
        J(2:4,2:4) = J(2:4,2:4) + J2;
        J(3:5,3:5) = J(3:5,3:5) + J3;
    elseif K == 4
        J = sym(zeros(6,6));

        J(1:3,1:3) = J(1:3,1:3) + J1;
        J(2:4,2:4) = J(2:4,2:4) + J2;
        J(3:5,3:5) = J(3:5,3:5) + J3;
        J(4:6,4:6) = J(4:6,4:6) + J4;
    elseif K == 5
        J = sym(zeros(7,7));

        J(1:3,1:3) = J(1:3,1:3) + J1;
        J(2:4,2:4) = J(2:4,2:4) + J2;
        J(3:5,3:5) = J(3:5,3:5) + J3;
        J(4:6,4:6) = J(4:6,4:6) + J4;
        J(5:7,5:7) = J(5:7,5:7) + J5;
    else
        J = sym(zeros(8,8));

        J(1:3,1:3) = J(1:3,1:3) + J1;
        J(2:4,2:4) = J(2:4,2:4) + J2;
        J(3:5,3:5) = J(3:5,3:5) + J3;
        J(4:6,4:6) = J(4:6,4:6) + J4;
        J(5:7,5:7) = J(5:7,5:7) + J5;
        J(6:8,6:8) = J(6:8,6:8) + J6;

    end
  

    %% Levi-Civita matrix

    Nlc = 3;
    [mats{1:Nlc}] = ndgrid(1:Nlc);
    pairsIndex = nchoosek(1:Nlc,2);
    lcMat = sign(prod(cat(Nlc+1,mats{pairsIndex(:,2)})-...
        cat(Nlc+1,mats{pairsIndex(:,1)}),Nlc+1));

    lcMat = sym(lcMat);

    %% Jacobi identity

    N = K + 2;

    Jacobiijkm = 0;

    for i = 1:N
        for j = 1:N
            for k = 1:N

                % alternating tensor
                ijk = [i, j, k];
                uijk = unique(ijk);

                if numel(uijk) <= 2
                    eijk = 0;
                else
                    [sijk,isort] = sort(ijk,'ascend');
                    eijk = lcMat(isort(1),isort(2),isort(3));
                end

                for m = 1:N
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
                    else
                        dJ_jk_xm = diff(J_jk,x8);
                    end

                    Jacobiijkm = Jacobiijkm + eijk*J_im*dJ_jk_xm;

                end
            end
        end
    end

    Jacobiijkmlist(K) = simplify(Jacobiijkm);  

end


%% Matrix of additional conditions for Hamiltonian structure, as new gyrostat is added

newcondmat = sym(NaN(6,1));

newcondmat(1,1) = simplify(Jacobiijkmlist(1));
for i = 2:6
    newcondmat(i,1) = simplify(Jacobiijkmlist(i)-Jacobiijkmlist(i-1));
end

