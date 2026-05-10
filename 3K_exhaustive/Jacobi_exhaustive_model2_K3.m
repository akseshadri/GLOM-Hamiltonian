%% Jacobi_exhaustive_model2_K3.m
%%
%% PURPOSE
%%   Exhaustive symbolic computation of the Jacobi condition for all
%%   3^3 = 27 representations of J for the K=3 sparse hierarchy
%%   (gyrostat 1 on {x1,x2,x3}, gyrostat 2 on {x3,x4,x5}, gyrostat 3
%%   on {x5,x6,x7}, M=7).
%%
%% PAPER REFERENCE
%%   Seshadri & Lakshmivarahan (2026a)
%%
%% METHOD
%%   Each gyrostat admits three representations L_A, L_B, L_C for its
%%   contribution to J, placing the state-dependent entries in the 3rd,
%%   2nd, or 1st column/row of its 3x3 block - respectively. 
%%
%% REQUIREMENTS: MATLAB Symbolic Math Toolbox
%%

clear, close all

p1 = sym('p1'); q1 = sym('q1'); r1 = sym('r1'); a1 = sym('a1'); b1 = sym('b1'); c1 = sym('c1');
p2 = sym('p2'); q2 = sym('q2'); r2 = sym('r2'); a2 = sym('a2'); b2 = sym('b2'); c2 = sym('c2');
p3 = sym('p3'); q3 = sym('q3'); r3 = sym('r3'); a3 = sym('a3'); b3 = sym('b3'); c3 = sym('c3');

%% State variables
x1 = sym('x1'); x2 = sym('x2'); x3 = sym('x3'); x4 = sym('x4'); x5 = sym('x5'); x6 = sym('x6'); x7 = sym('x7');

%% Define each gyrostat's contribution to J

% gyrostat 1
J11 = [0 -c1 p1*x2+b1; ...  % x1, x2, x3
    c1 0 q1*x1-a1; ...
    -(p1*x2+b1) -(q1*x1-a1) 0];

J12 = [0 p1*x3-c1 b1; ...
    -p1*x3+c1 0 -r1*x1-a1; ...
    -b1 r1*x1+a1 0];

J13 = [0 -q1*x3-c1 -r1*x2+b1; ...
    q1*x3+c1 0 -a1 ; ...
    r1*x2-b1 a1 0];

% gyrostat 2
J21 = [0 -c2 p2*x4+b2; ...  % x3, x4, x5
    c2 0 q2*x3-a2; ...
    -(p2*x4+b2) -(q2*x3-a2) 0];

J22 = [0 p2*x5-c2 b2; ...
    -p2*x5+c2 0 -r2*x3-a2; ...
    -b2 r2*x3+a2 0];

J23 = [0 -q2*x5-c2 -r2*x4+b2; ...
    q2*x5+c2 0 -a2 ; ...
    r2*x4-b2 a2 0];

% gyrostat 3
J31 = [0 -c3 p3*x6+b3; ...  % x5, x6, x7
    c3 0 q3*x5-a3; ...
    -(p3*x6+b3) -(q3*x5-a3) 0];

J32 = [0 p3*x7-c3 b3; ...
    -p3*x7+c3 0 -r3*x5-a3; ...
    -b3 r3*x5+a3 0];

J33 = [0 -q3*x7-c3 -r3*x6+b3; ...
    q3*x7+c3 0 -a3; ...
    r3*x6-b3 a3 0];


J1all = cat(3,J11,J12,J13);
J2all = cat(3,J21,J22,J23);
J3all = cat(3,J31,J32,J33);


%% Check Jacobi condition

count = 0;

Jacobiijkmlist = sym(zeros(3^3,1));

for d1 = 1:3
    for d2 = 1:3

        for d3 = 1:3
            count = count + 1

            if count == 1
                dlist = [d1,d2,d3];
            else
                dlist = cat(1,dlist,[d1 d2 d3]);

            end

            J1 = squeeze(J1all(:,:,d1));
            J2 = squeeze(J2all(:,:,d2));
            J3 = squeeze(J3all(:,:,d3));


            K = 3;


            J = sym(zeros(7,7));

            i1 = [1 2 3];
            J(i1,i1) = J(i1,i1) + J1;

            i2 = [3 4 5];
            J(i2,i2) = J(i2,i2) + J2;

            i3 = [5 6 7];
            J(i3,i3) = J(i3,i3) + J3;

            %% Levi-Civita matrix

            Nlc = 3;
            [mats{1:Nlc}] = ndgrid(1:Nlc);
            pairsIndex = nchoosek(1:Nlc,2);
            lcMat = sign(prod(cat(Nlc+1,mats{pairsIndex(:,2)})-...
                cat(Nlc+1,mats{pairsIndex(:,1)}),Nlc+1));

            lcMat = sym(lcMat);

            %% Jacobi identity

            N = 7;

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
                            elseif m == 8
                                dJ_jk_xm = diff(J_jk,x8);
                            else
                                dJ_jk_xm = diff(J_jk,x9);

                            end

                            Jacobiijkm = Jacobiijkm + eijk*J_im*dJ_jk_xm;

                        end
                    end
                end
            end

            Jacobiijkmlist(count) = simplify(Jacobiijkm);


        end
    end
end



