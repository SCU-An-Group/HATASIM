clear;
clc;
close all;
%% Step 1:
[Coords, Ebc, Loads, Connect, E, A, I] = model();

%% Step 2: 
[Model, Node, Elem, P] = pre_processing(Coords, Ebc, Loads, Connect, E, A, I);

%% Step 3:                                               
U   = zeros(Model.neq, 1);   
d_U = zeros(Model.neq, 1);   

%====================================================================================================================================
% Start incremental process
tol = 1e-8;       
maxIter = 100;     
res1 = []; 
res2 = [];  

iter = 0;

tic;

while iter < maxIter
    iter = iter + 1;
    fprintf('iteration %d\n', iter); 
    
    [Kt, Elem] = tangStiffMtx(Model, Elem, U);
    [F,  Elem] = intForces(Model, Elem);
    
    if checkSingularMtx(Model, Kt)
        fprintf(1, 'Singular tangent matrix!\nStep = %d\nIter = %d\n', iter);
        break;
    end
    
    R = P - F;
    res1(end+1) = abs(R(1));
    res2(end+1) = abs(R(2));
    erro = norm(R(1:Model.neqf)) / norm(P(1:Model.neqf));
    
    if erro < tol
        disp('Final displacement vector U:');
        disp(U);
        break;
    end
        
    d_U = solveLinearSystem(Model, Kt, R);
    U   = U + d_U;
end

elapsed_time = toc;  
fprintf('Total runtime: %.4f seconds\n', elapsed_time);

figure;
plot(1:length(res1), res1, 'b-o', 'LineWidth', 1.5); hold on;
if exist('res2', 'var')
    plot(1:length(res2), res2, 'r--s', 'LineWidth', 1.5);
    legend('|R(1)|','|R(2)|');
else
    legend('|R(1)|');
end
xlabel('Iteration');
ylabel('Residual');
title('Convergence history');
grid on;

if iter == maxIter
    fprintf('Warning: Newton iterations did not converge at load step %d.\n', iter);
end
    

origX = Coords(:,1);
origY = Coords(:,2);
defX = zeros(Model.nnp,1);
defY = zeros(Model.nnp,1);
for i = 1:Model.nnp
    dof = Node(i).dof;      
    defX(i) = Coords(i,1) + U(dof(1));
    defY(i) = Coords(i,2) + U(dof(2));
end

assignin('base', 'defX', defX);
assignin('base', 'defY', defY);
figure;
plot(origX, origY, 'k-o','LineWidth',1.5); hold on;
plot(defX, defY, 'r-o','LineWidth',1.5);
XCOOR_all = []; 
for i = 1:Model.nel
    if isfield(Elem(i), 'XCOOR') && ~isempty(Elem(i).XCOOR)
        plot(Elem(i).XCOOR(:,1), Elem(i).XCOOR(:,2), 'b-', 'LineWidth', 1.2);
        XCOOR_all = [XCOOR_all; Elem(i).XCOOR];  % 累加每个元素的XCOOR
    end
end

XCOOR_all = unique(XCOOR_all, 'rows');

assignin('base', 'XCOOR_all', XCOOR_all);

assignin('base', 'XCOOR_all', XCOOR_all);
legend('Original','Deformed');
xlabel('X'); ylabel('Y');
title('Structure Shape Comparison'); grid on;
 
allTensions = []; 
for i = 1:Model.nel
    if isfield(Elem(i), 'Tension') && ~isempty(Elem(i).Tension)
        allTensions = [allTensions; Elem(i).Tension(:)];
    end
end

allTensions = unique(allTensions, 'stable');

assignin('base', 'allTensions', allTensions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Coords,Ebc,Loads,Connect,E,A,I] = model

    L = 999.9;        
    nElem = 10;     
    nNode = nElem+1; 

    Coords = zeros(nNode, 2);
    for i = 0:nElem
        Coords(i+1, :) = [0.0, L*i/nElem];
    end
    
    Ebc = repmat([0 0 1], nNode, 1);  
    Ebc(1, :) = [1 1 1];              


    Loads = [
      [0.99  -0.99454  0.0];
      [2.96901  -2.73708  0.0];
      [4.8883  -3.97823  0.0];
      [6.78907  -4.83284  0.0];
      [8.75855  -5.43688  0.0];
      [10.77131  -5.80626  0.0];
      [12.787  -5.94787  0.0];
      [14.76398  -5.87206  0.0];
      [16.66412  -5.5948  0.0];
      [18.45533  -5.13712  0.0];
      [9.66+118.4608  -2.433-197.876+722.6637  0.0];
    ];
    
    Connect = [(1:nElem)'  (2:nElem+1)'];
  
    
    E(1:nElem)  = 95e9;
    A(1:nElem)  = 2.82743E-05;
    I(1:nElem)  = 3.641573376e-3;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Model,Node,Elem, P] = pre_processing(Coords,Ebc,Loads,Connect,E,A,I)
    nnp = size(Coords,1);
    nel = size(Connect,1);
    neq = 3 * nnp;

    ID = zeros(nnp,3);
    neqfixed = 0;
    for i = 1:nnp
        for j = 1:3
            if (Ebc(i,j) == 1)
                neqfixed = neqfixed + 1;
                ID(i,j) = 1;
            end
        end
    end

    neqfree = neq - neqfixed;
    countS = neqfree;
    countF = 0;
    for i = 1:nnp
        for j = 1:3
            if ID(i,j) == 0
                countF = countF + 1;
                ID(i,j) = countF;
            else
                countS = countS + 1;
                ID(i,j) = countS;
            end
        end
    end

    % Assemble elements gather vectors
    GLE = zeros(nel, 6);
    for i = 1:nel
        for j = 1:6
            if (j <= 3)
                GLE(i,j) = ID(Connect(i,1),j);
            else
                GLE(i,j) = ID(Connect(i,2),j-3);
            end
        end
    end

    disp('GLE matrix:');
    disp(GLE);
    disp('ID matrix:');
    disp(ID);

    % Assemble load vector
    P  = zeros(neq,1);
    for i = 1:nnp
        for j = 1:3
            P(ID(i,j)) = P(ID(i,j)) + Loads(i,j);
        end
    end

    % Create model structure
    Model = struct('nnp',nnp,'nel',nel,'neq',neq,'neqf',neqfree,'neqc',neqfixed);

    % Create vector of node structures
    Node(Model.nnp,1) = struct('x',[],'y',[],'px',[],'py',[],'mz',[],'dof',[]);
    for i = 1:Model.nnp
        Node(i).x   = Coords(i,1);
        Node(i).y   = Coords(i,2);
        Node(i).px  = Loads(i,1);
        Node(i).py  = Loads(i,2);
        Node(i).mz  = Loads(i,3);
        Node(i).dof = ID(i,:);
    end

    % Create vector of element structures
    Elem(Model.nel,1) = struct('n1',[],'n2',[],'E',[],'A',[],'I',[],...
                               'fi',[],'fi_1',[],'fn',[],'fn_1',[],...
                               'gle',[],'ke',[],'kt',[]);
    for i = 1:Model.nel
        % Get element node numbers
        n1 = Node(Connect(i,1));
        n2 = Node(Connect(i,2));
    
        % Set element properties
        Elem(i).n1      = n1;
        Elem(i).n2      = n2;
        Elem(i).E       = E(i);
        Elem(i).A       = A(i);
        Elem(i).I       = I(i);
        Elem(i).fi      = zeros(6,1);
        Elem(i).fi_1    = zeros(6,1);
        Elem(i).fn      = zeros(3,1);
        Elem(i).fn_1    = zeros(3,1);
        Elem(i).gle     = GLE(i,:);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Kt,Elem] = tangStiffMtx(Model,Elem,U)
    Kt = zeros(Model.neq,Model.neq);

    for i = 1:Model.nel         
        [ke, F_loc, XCOOR, TensionVec] = cableStiffMtx(Elem(i), U);
            
        Elem(i).fi = F_loc;
        Elem(i).XCOOR = XCOOR;
        Elem(i).Tension = TensionVec;

        kg = zeros(6,6); 
        kt = ke + kg ;
        
        Elem(i).ke = ke;

        k = kt;

        gle = Elem(i).gle;
        Kt(gle,gle) = Kt(gle,gle) + k;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function D = solveLinearSystem(Model,K,P)
    Kff = K(1:Model.neqf, 1:Model.neqf);
    Pf  = P(1:Model.neqf);
    Ds  = zeros(Model.neqc,1);
    Df = Kff \ Pf;
    D = [Df; Ds];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F,Elem] = intForces(Model,Elem)
    F = zeros(Model.neq,1);
    for i = 1:Model.nel

        fg = Elem(i).fi;
        
        gle = Elem(i).gle;
        F(gle) = F(gle) + fg;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function singular = checkSingularMtx(Model,K)
    singular = 0;
    if (rcond(K(1:Model.neqf,1:Model.neqf)) < 10e-12)
        singular = 1;
    end
end
