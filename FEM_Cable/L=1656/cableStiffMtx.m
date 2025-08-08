function [ke,F_loc, XCOOR, TensionVec] = cableStiffMtx(Elem, U)

    E = Elem.E;
    A = Elem.A;

    Node1 = Elem.n1;
    Node2 = Elem.n2;

    % Nodal displacements
    dx1 = U(Node1.dof(1));
    dy1 = U(Node1.dof(2));
    dx2 = U(Node2.dof(1));
    dy2 = U(Node2.dof(2));

    % New nodal coordinates
    x1 = Node1.x + dx1;
    y1 = Node1.y + dy1;
    x2 = Node2.x + dx2;
    y2 = Node2.y + dy2;


    Z = x2 - x1;
    T = y2 - y1;

    MEMNO = 1;
    TEMP = 0.0;       
    ET = 0.65e-5;       
    XL0 =1656.7/10;
    W0 = 0.269;          
    IMPRE = 0;          
    NPTS = 40;          
    a = 1e-4;           
    
    [ke,F_loc, XCOOR, TensionVec] = CableTSM(MEMNO, Z, T, A, E, TEMP, ET, XL0, W0, IMPRE, NPTS, a);

    P1_global = [x1; y1];
    XCOOR = XCOOR + P1_global(:)';

end