function OUT  = weighted_average_outputs(OUT_A, OUT_B, weight_A, weight_B)
    % Calculate weighted averages for each field
    OUT.Ea = weight_A * OUT_A.Ea + weight_B * OUT_B.Ea;
    OUT.QF = weight_A * OUT_A.QF + weight_B * OUT_B.QF;
    OUT.R = weight_A * OUT_A.R + weight_B * OUT_B.R;
    OUT.QS = weight_A * OUT_A.QS + weight_B * OUT_B.QS;
    OUT.QT = weight_A * OUT_A.QT + weight_B * OUT_B.QT;
    OUT.Sf = weight_A * OUT_A.Sf + weight_B * OUT_B.Sf;
    OUT.Su = weight_A * OUT_A.Su + weight_B * OUT_B.Su;
    OUT.Ss = weight_A * OUT_A.Ss + weight_B * OUT_B.Ss;
    OUT.St = weight_A * OUT_A.St + weight_B * OUT_B.St;
    OUT.AL = weight_A * OUT_A.AL + weight_B * OUT_B.AL;
    OUT.IE = weight_A * OUT_A.IE + weight_B * OUT_B.IE;
    OUT.SE = weight_A * OUT_A.SE + weight_B * OUT_B.SE;
    OUT.Ei = weight_A * OUT_A.Ei + weight_B * OUT_B.Ei;
    OUT.Et = weight_A * OUT_A.Et + weight_B * OUT_B.Et;
    OUT.S_canopy = weight_A * OUT_A.S_canopy + weight_B * OUT_B.S_canopy;
    OUT.pot_inf = weight_A * OUT_A.pot_inf + weight_B * OUT_B.pot_inf;
    OUT.P        = OUT_A.P;
    OUT.Ep        = OUT_A.Ep;

end