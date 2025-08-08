function [ke, F_loc, XCOOR, TensionVec] = CableTSM(MEMNO, Z, T, A, E, TEMP, ET, XLO, WO, IMPRE, NPTS, a)


[F03, F04, H, V, F_loc, XCOOR, TensionVec] = pcaFX2(MEMNO, Z, T, A, E, TEMP, ET, XLO, WO, IMPRE, NPTS);


[F03_pp, F04_pp] = pcaFX2(MEMNO, Z, T*(1+a), A, E, TEMP, ET, XLO, WO, IMPRE, NPTS);

aaa = 1e-3;
delta = aaa * XLO;  


if abs(T) < 1e-10
    [F03_pp, F04_pp]       = pcaFX2(MEMNO, Z, T+delta,  A, E, TEMP, ET, XLO, WO, IMPRE, NPTS);

    H = Z;
    V = T;

    a3_val = (F03_pp - F03)/delta;
    a4_val = (F04_pp - F04)/delta;

else
    [F03_pp, F04_pp]       = pcaFX2(MEMNO, Z, T*(1+a),  A, E, TEMP, ET, XLO, WO, IMPRE, NPTS);


    H = Z;
    V = T;

    a3_val = (F03_pp - F03)/(a * V);
    a4_val = (F04_pp - F04)/(a * V);
end




if abs(Z) < 1e-10
    [F03_prime, F04_prime]  = pcaFX2(MEMNO, Z+delta, T,  A, E, TEMP, ET, XLO, WO, IMPRE, NPTS);

    H = Z;
    V = T;


    a1_val = (F03_prime - F03) / (delta);
    a2_val = (F04_prime - F04) / (delta);

else
    [F03_prime, F04_prime]       = pcaFX2(MEMNO, Z*(1+a), T,  A, E, TEMP, ET, XLO, WO, IMPRE, NPTS);


    H = Z;
    V = T;

    a1_val = (F03_prime - F03) / (a * H);
    a2_val = (F04_prime - F04) / (a * H);
end


ke = [  a1_val,  a3_val,  0,      -a1_val, -a3_val,  0;
       a2_val,  a4_val,  0,      -a2_val, -a4_val,  0;
       0,       0,       0,       0,       0,       0;
      -a1_val, -a3_val,  0,       a1_val,  a3_val,   0;
      -a2_val, -a4_val,  0,       a2_val,  a4_val,   0;
       0,       0,       0,       0,       0,       0];

end