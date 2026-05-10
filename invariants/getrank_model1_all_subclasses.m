%% getrank_model1_all_subclasses.m
%%
%% PURPOSE
%%   Enumerate all 2^12 parameter subclasses of Model 1 (K=2, M=4,
%%   dense coupling) and, for each, compute: the invariant count,
%%   whether the subclass is Hamiltonian, the rank of J, and the
%%   explicit Casimir if Hamiltonian. 

%% PAPER REFERENCE
%%   Seshadri & Lakshmivarahan (2026a)
%%
%% METHOD
%%   Each of the 12 binary flags {ip1,iq1,...,ic2} indicates whether
%%   the corresponding parameter is included (=1) or set to zero (=0).
%%   For each valid subclass (satisfying energy-conservation constraints),
%%   the constraint matrix A is built symbolically, its rank computed, and the Jacobi condition evaluated
%%   for the canonical J matrix.
%%
%% REQUIREMENTS: MATLAB Symbolic Math Toolbox, Statistics Toolbox
%%
%% KEY OUTPUTS
%%   count_inv   -- invariant count per subclass
%%   isHam       -- Hamiltonian flag (0/1) per subclass
%%   rankJ       -- rank of J per subclass
%%   Casimir     -- symbolic Casimir expression (Hamiltonian cases only)
%%

%% Cases with at most one gyrostat having 2 nonlinear terms
%% SI Fig 1 and 2 are plotted here


clear, close all hidden

%% Loop across subclasses

count = 0;

for i = 0:2^12-1

    % Parameters
    p1 = sym('p1'); q1 = sym('q1'); r1 = sym('r1');
    p2 = sym('p2'); q2 = sym('q2'); r2 = sym('r2');

    a1 = sym('a1'); b1 = sym('b1'); c1 = sym('c1');
    a2 = sym('a2'); b2 = sym('b2'); c2 = sym('c2');

    % Variables
    x1 = sym('x1'); x2 = sym('x2'); x3 = sym('x3'); x4 = sym('x4');

    isincluded = dec2bin(i,12);

    % assign inclusion of terms
    ip1 = isequal(isincluded(1),'1')*1;
    iq1 = isequal(isincluded(2),'1')*1;
    ir1 = isequal(isincluded(3),'1')*1;
    ia1 = isequal(isincluded(4),'1')*1;
    ib1 = isequal(isincluded(5),'1')*1;
    ic1 = isequal(isincluded(6),'1')*1;
    ip2 = isequal(isincluded(7),'1')*1;
    iq2 = isequal(isincluded(8),'1')*1;
    ir2 = isequal(isincluded(9),'1')*1;
    ia2 = isequal(isincluded(10),'1')*1;
    ib2 = isequal(isincluded(11),'1')*1;
    ic2 = isequal(isincluded(12),'1')*1;

    nnonlin1 = ip1+iq1+ir1;
    nnonlin2 = ip2+iq2+ir2;

    nlin = ia1+ib1+ic1+ia2+ib2+ic2;


   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
   % Determine if legitimate model

   if ((nnonlin1 == 2) & (nnonlin2 == 3)) | ((nnonlin1 == 3) & (nnonlin2 == 2)) | ((nnonlin1 == 3) & (nnonlin2 == 3))
        count = count + 1

        % Save subclass
        pilistc = [ip1 iq1 ir1 ia1 ib1 ic1 ip2 iq2 ir2 ia2 ib2 ic2];

        if count == 1
            pilist = pilistc;
        else
            pilist = cat(1,pilist,pilistc);
        end

        % Assign parameters

        if nnonlin1 == 3

            p1 = ip1*p1;
            q1 = iq1*q1;
            r1 = ir1*(-1)*(p1 + q1);

        else % two nonlinear terms

            if ir1 == 0
                p1 = ip1*p1;
                q1 = -p1;
                r1 = 0;

            elseif iq1 == 0
                p1 = ip1*p1;
                q1 = 0;
                r1 = -p1;

            elseif ip1 == 0
                p1 = 0;
                q1 = iq1*q1;
                r1 = -q1;
            end
        end

        if nnonlin2 == 3

            p2 = ip2*p2;
            q2 = iq2*q2;
            r2 = ir2*(-1)*(p2 + q2);

        else % two nonlinear terms

            if ir2 == 0
                p2 = ip2*p2;
                q2 = -p2;
                r2 = 0;

            elseif iq2 == 0
                p2 = ip2*p2;
                q2 = 0;
                r2 = -p2;

            elseif ip2 == 0
                p2 = 0;
                q2 = iq2*q2;
                r2 = -q2;
            end
        end

        a1 = ia1*a1;
        b1 = ib1*b1;
        c1 = ic1*c1;

        a2 = ia2*a2;
        b2 = ib2*b2;
        c2 = ic2*c2;

        % Matrix J
        J1 = [0 -c1 p1*x2+b1; ...
            c1 0 q1*x1-a1; ...
            -(p1*x2+b1) -(q1*x1-a1) 0];

        J2 = [0 -c2 p2*x3+b2; ...
            c2 0 q2*x2-a2; ...
            -(p2*x3+b2) -(q2*x2-a2) 0];

        J = sym(zeros(4,4));

        J(1:3,1:3) = J(1:3,1:3) + J1;
        J(2:4,2:4) = J(2:4,2:4) + J2;

        clear J1 J2 

        % Full matrix
        A = [p1 q1 r1 0 r2 0 0 0 0 ; ...
            -c1 c1 0 0 -b2 0 0 r1 0 ; ...
            0 -(a1+c2) (a1+c2) 0 0 p1 0 0 r2; ...
            b1 0 -b1 0 a2 0 q1 0 0 ; ...
            0 p2 q2 r2 p1 0 0 0 0 ; ...
            0 0 -a2 a2 b1 0 p2 0 0; ...
            0 b2 0 -b2 -c1 0 0 q2 0; ...
            0 0 0 0 0 0 c1 -b1 0; ...
            0 0 0 0 0 -c1 0 (a1+c2) -b2; ...
            0 0 0 0 0 b1 -(a1+c2) 0 a2; ...
            0 0 0 0 0 0 b2 -a2 0];
    

        % Matrix rank
        R_A = rank(A);

        if count == 1
            rAlist = R_A;
            Astack = A;
        else
            rAlist = cat(1,rAlist,R_A);
            Astack = cat(3,Astack,A);
        end

        if count == 1
            Jstack = J;
        else
            Jstack = cat(3,Jstack,J);
        end

        % Store list of parameters
        parami = [p1 q1 r1 a1 b1 c1 p2 q2 r2 a2 b2 c2];
        if count == 1
            parammat = parami;
        else
            parammat = cat(1,parammat,parami);
        end

        % Store binary representation
        if count == 1
            binrep = isincluded;
        else
            binrep = cat(1,binrep,isincluded);
        end

        % Store array representation
        arrayparami = [ip1 iq1 ir1 ia1 ib1 ic1 ip2 iq2 ir2 ia2 ib2 ic2];
        if count == 1
            binrepn = arrayparami;
        else
            binrepn = cat(1,binrepn,arrayparami);
        end       

        % Is the subclass Hamiltonian?
        Htest1 = p1*p2;
        Htest2 = p2*b1 - p1*b2 - q2*c1;

        if isequal(Htest1,sym(0)) && isequal(Htest2,sym(0)) % Hamiltonian



            if count == 1
                Hlist = 1;
            else
                Hlist = cat(1,Hlist,1);
            end



        else
            if count == 1
                Hlist = 0;
            else
                Hlist = cat(1,Hlist,0);
            end
        end

        % Store Casimir
        rJi = rank(J);

        % Nullspace
        NJi = null(J); dimNJ = 4-rJi;

        NJmati = sym(zeros(4,4));
        if dimNJ > 0
            NJmati(:,1:dimNJ) = NJi;
        end

        if count == 1
            RJlist = rJi;

            NJmat = NJmati;

        else
            RJlist = cat(1,RJlist,rJi);
            NJmat = cat(3,NJmat,NJmati);
        end


   end

end


%% Regression tree for Hamiltonian structure
y1 = Hlist;

% fit tree
qinvtree = fitrtree(binrepn,y1,'MinParentSize',5,'PredictorNames',{'p1','q1','r1','a1','b1','c1','p2','q2','r2','a2','b2','c2'});

% view tree
view(qinvtree,'Mode','graph')

%% Regression tree for number of invariants

y2 = NaN(numel(rAlist),1);

% prepare pedictor
i5 = find(rAlist == 5);
y2(i5) = 3;

i7 = find(rAlist == 7);
y2(i7) = 2;

i8 = find(rAlist == 8);
y2(i8) = 1;

qinvtree = fitrtree(binrepn,y2,'MinParentSize',5,'PredictorNames',{'p1','q1','r1','a1','b1','c1','p2','q2','r2','a2','b2','c2'});

% view tree
view(qinvtree,'Mode','graph')

%% Regression tree for number of invariants among Hamiltonian models
iH = find(Hlist == 1);

y3 = y2(iH);
binrepnH = binrepn(iH,:);

qinvtree = fitrtree(binrepnH,y3,'MinParentSize',5,'PredictorNames',{'p1','q1','r1','a1','b1','c1','p2','q2','r2','a2','b2','c2'});

% view tree
view(qinvtree,'Mode','graph')


%% Analysis

figure, 
subplot(3,1,1), plot(y2,'k+'), grid on, ylabel('# invariants')
subplot(3,1,2), plot(y1,'k+'), grid on, ylabel('Hamiltonian')
subplot(3,1,3), plot(RJlist,'k+'), grid on, ylabel('rank(J)')

figure,
subplot(1,2,1), plot(y1,y2,'ro','MarkerSize',10), xlabel('Hamiltonian'), ylabel('# invariants'), grid on
subplot(1,2,2), plot(y2(iH),RJlist(iH),'ro','MarkerSize',10), xlabel('# invariants (Hamiltonian)'), ylabel('rank(J)'), grid on

%% Null space for nondegenerate Hamiltonian case, p2 = c1 = b2 = 0

strhamilt1 = [1 1 1 1 1 0 0 1 1 1 0 1];

for i = 1:size(binrepn,1)
    binrepni = binrepn(i,:);
    if binrepni == strhamilt1
        out1 = i;
    end
end

Aout1 = squeeze(Astack(:,:,out1));
Jout1 = squeeze(Jstack(:,:,out1));

InvAout1 = null(Aout1);
NJout1 = null(Jout1);

%% Null space for degenerate Hamiltonian case, p1 = b1 = c1 = 0

strhamilt2 = [0 1 1 1 0 0 1 1 1 1 1 1];

for i = 1:size(binrepn,1)
    binrepni = binrepn(i,:);
    if binrepni == strhamilt2
        out2 = i;
    end
end

Aout2 = squeeze(Astack(:,:,out2));
Jout2 = squeeze(Jstack(:,:,out2));

InvAout2 = null(Aout2);
NJout2 = null(Jout2);

%% Null space for degenerate Hamiltonian case, p2 = c1 = b2 = 0 and b1 = a2 = 0

strhamilt3 = [1 1 1 1 0 0 0 1 1 0 0 1];

for i = 1:size(binrepn,1)
    binrepni = binrepn(i,:);
    if binrepni == strhamilt3
        out3 = i;
    end
end

Aout3 = squeeze(Astack(:,:,out3));
Jout3 = squeeze(Jstack(:,:,out3));

InvAout3 = null(Aout3);
NJout3 = null(Jout3);

%% Null space for all Hamiltonian cases
iC2 = find(y2 == 2 & RJlist == 2 & Hlist == 1);
iC3 = find(y2 == 3 & RJlist == 2 & Hlist == 1);


% 2 invariants
for i = iC2'
    Ai = squeeze(Astack(:,:,i));
    Ji = squeeze(Jstack(:,:,i));

    InvAi = null(Ai)
    NJi = null(Ji)
end

% 3 invariants
for i = iC3'
    Ai = squeeze(Astack(:,:,i));
    Ji = squeeze(Jstack(:,:,i));

    InvAi = null(Ai)
    NJi = null(Ji)
end