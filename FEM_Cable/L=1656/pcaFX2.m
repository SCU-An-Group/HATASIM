function [F03, F04, H, V, F_loc, XCOOR, TensionVec] = pcaFX2(MEMNO, Z, T, A, E, TEMP, ET, XLO, WO, IMPRE, NPTS)

EPS = 1e-6;
PI  = pi;



H = Z;
V = T;
ICODE = 3; 


XL = XLO * (1 + ET * TEMP);
W = WO * XLO / XL;


LL = 0;

KK = 0;
if V > 0
    KK = 1;
    V = -V;
    H = -H;
end

CORD = sqrt(H^2 + V^2);
AMBDA = 1e6;
if H ~= 0
    if XL > CORD
        AMBDA = sqrt( ((XL^2 - V^2) / (H^2) - 1) * 3 );
    else
        AMBDA = 0.20;
    end
end

F01 = -(W * H) / (2 * AMBDA);
COT = 1 / tanh(AMBDA);
F02 = (W/2) * (-V * COT + XL);
DF1 = 0;
DF2 = 0;

maxIter = 1000;
iter = 0;
while true
    F01 = F01 + DF1;
    F02 = F02 + DF2;
    
    F04 = W * XL - F02;
    F03 = -F01;
    
    TI = sqrt(F01^2 + F02^2);
    TJ = sqrt(F03^2 + F04^2);
    
    F_val = F04 + TJ;
    FF = TI - F02;
    if FF < 1e-4
        FF = 1e-4;
    end
    G = F_val / FF;
    if G < 1e-4
        G = 1e-4;
    end
    
    AAH = (1/W)*log(G) + XL/(A*E);
    AH = -F01 * AAH;
    BV = (TJ^2 - TI^2)/(2*E*A*W) + (TJ - TI)/W;
    
    CA = H - AH;
    CB = V - BV;
    if (abs(CA) <= EPS) && (abs(CB) <= EPS)
        XLAFST = XL + (F04*TJ + F02*TI + F01^2 * log(G))/(2*E*A*W);
        break;
    end
    
    if TJ < 1e-4
        TJ = 1e-4;
    end
    
    B2 = -(F02/TI + F04/TJ)/W - XL/(E*A);
    A1 = -AAH - B2 - XL/(E*A);
    A2 = (F01/W)*(1/TJ - 1/TI);
    DET = A1 * B2 - A2^2;
    DF1 = (CA * B2 - CB * A2) / DET;
    DF2 = (A1 * CB - A2 * CA) / DET;
    
    iter = iter + 1;
    if iter > maxIter
        warning('pcaFX2 未能在最大迭代次数内收敛！');
        break;
    end
end

FOC = zeros(6,1);
FOC(1) = F01 * (1 - 2 * KK);
FOC(2) = F02 + KK * (F04 - F02);
FOC(4) = F03 * (1 - 2 * KK);
FOC(5) = F04 + KK * (F02 - F04);

F03 = FOC(4);
F04 = FOC(5);


%
SUBXL = XLO / (NPTS - 1);
XL_tmp = -SUBXL;
TensionVec = zeros(NPTS+1,1);

for MM = 1:NPTS
    XL_tmp = XL_tmp + SUBXL;
    F04_p = W * XL_tmp - F02;
    F03_p = -F01;
    TI_p = sqrt(F01^2 + F02^2);
    TJ_p = sqrt(F03_p^2 + F04_p^2);


    if MM == 1
        TensionVec(1) = TI_p;
    end

    TensionVec(MM+1) = TJ_p;

    F_p = F04_p + TJ_p;
    FF_p = TI_p - F02;
    if FF_p < 1e-4
        FF_p = 1e-4;
    end
    G_p = F_p / FF_p;
    if G_p < 1e-4
        G_p = 1e-4;
    end
    AAH_p = (1 / W) * log(G_p) + XL_tmp / (A * E);
    AH_p = -F01 * AAH_p;
    BV_p = (TJ_p^2 - TI_p^2) / (2 * E * A * W) + (TJ_p - TI_p) / W;

    MN = MM + (NPTS - 2 * MM + 1) * KK;
    XCOOR(MN, 1) = AH_p + Z * KK;
    XCOOR(MN, 2) = BV_p + T * KK;
end

TensionVec = flipud(TensionVec);
F_loc = [FOC(1);FOC(2);0;FOC(4);FOC(5);0];

end
