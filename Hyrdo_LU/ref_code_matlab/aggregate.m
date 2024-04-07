function OUT = aggregate(OUT, dates)

% Aggregate Monthly
[E_mon, ~]       = make_monthly_new(OUT.Ea, dates, 2);
[Ss_mon, ~]      = make_monthly_new(OUT.Ss, dates, 1);
[Su_mon, ~]      = make_monthly_new(OUT.Su, dates, 1);
[St_mon, ~]      = make_monthly_new(OUT.St, dates, 1);
[QS_mon, ~]      = make_monthly_new(OUT.QS, dates, 2);
[QF_mon, ~]      = make_monthly_new(OUT.QF, dates, 2);
[QT_mon, ~]      = make_monthly_new(OUT.QT, dates, 2);
[IE_mon, ~]      = make_monthly_new(OUT.IE, dates, 2);
[SE_mon, ~]      = make_monthly_new(OUT.SE, dates, 2);
[R_mon, ~]       = make_monthly_new(OUT.R, dates, 2);
[P_mon, ~]       = make_monthly_new(OUT.P, dates, 2);
[PET_mon, ~]     = make_monthly_new(OUT.Ep, dates, 2);
[Ei_mon, ~]      = make_monthly_new(OUT.Ei, dates, 2);
[Et_mon, ~]      = make_monthly_new(OUT.Et, dates, 2);
[S_canopy_mon, ~]= make_monthly_new(OUT.S_canopy, dates, 1);
[pot_inf_mon, ~] = make_monthly_new(OUT.pot_inf, dates, 2);

% Aggregate mean_monthly and mean annual
[E_mean_mon, E_mean_annual]         = make_means_new(OUT.Ea, dates);
[Ss_mean_mon, Ss_mean_annual]       = make_means(OUT.Ss, dates);
[Su_mean_mon, Su_mean_annual]       = make_means(OUT.Su, dates);
[St_mean_mon, St_mean_annual]       = make_means(OUT.St, dates);
[QS_mean_mon, QS_mean_annual]       = make_means_new(OUT.QS, dates);
[QF_mean_mon, QF_mean_annual]       = make_means_new(OUT.QF, dates);
[QT_mean_mon, QT_mean_annual]       = make_means_new(OUT.QT, dates);
[IE_mean_mon, IE_mean_annual]       = make_means_new(OUT.IE, dates);
[SE_mean_mon, SE_mean_annual]       = make_means_new(OUT.SE, dates);
[R_mean_mon, R_mean_annual]         = make_means_new(OUT.R, dates);
[P_mean_mon, P_mean_annual]         = make_means_new(OUT.P, dates);
[PET_mean_mon, PET_mean_annual]     = make_means_new(OUT.Ep, dates);
[Ei_mean_mon, Ei_mean_annual]       = make_means_new(OUT.Ei, dates);
[Et_mean_mon, Et_mean_annual]       = make_means_new(OUT.Et, dates);
[S_canopy_mean_mon, S_canopy_mean_annual] = make_means_new(OUT.S_canopy, dates);
[pot_inf_mean_mon, pot_inf_mean_annual] = make_means_new(OUT.pot_inf, dates);

% Save Monthly
OUT.E_mon           = E_mon;
OUT.Ss_mon          = Ss_mon;
OUT.Su_mon          = Su_mon;
OUT.St_mon          = St_mon;
OUT.QS_mon          = QS_mon;
OUT.QF_mon          = QF_mon;
OUT.QT_mon          = QT_mon;
OUT.IE_mon          = IE_mon;
OUT.SE_mon          = SE_mon;
OUT.R_mon           = R_mon;
OUT.P_mon           = P_mon;
OUT.PET_mon         = PET_mon;
OUT.Ei_mon          = Ei_mon;
OUT.Et_mon          = Et_mon;
OUT.S_canopy_mon    = S_canopy_mon;
OUT.pot_inf_mon     = pot_inf_mon;

% Save Mean Monthly
OUT.E_mean_mon          = E_mean_mon;
OUT.Ss_mean_mon         = Ss_mean_mon;
OUT.Su_mean_mon         = Su_mean_mon;
OUT.St_mean_mon         = St_mean_mon;
OUT.QS_mean_mon         = QS_mean_mon;
OUT.QF_mean_mon         = QF_mean_mon;
OUT.QT_mean_mon         = QT_mean_mon;
OUT.IE_mean_mon         = IE_mean_mon;
OUT.SE_mean_mon         = SE_mean_mon;
OUT.R_mean_mon          = R_mean_mon;
OUT.P_mean_mon          = P_mean_mon;
OUT.PET_mean_mon        = PET_mean_mon;
OUT.Ei_mean_mon         = Ei_mean_mon;
OUT.Et_mean_mon         = Et_mean_mon;
OUT.S_canopy_mean_mon   = S_canopy_mean_mon;
OUT.pot_inf_mean_mon    = pot_inf_mean_mon;

% Save Mean Annual
OUT.E_mean_annual       = E_mean_annual;
OUT.Ss_mean_annual      = Ss_mean_annual;
OUT.Su_mean_annual      = Su_mean_annual;
OUT.St_mean_annual      = St_mean_annual;
OUT.QS_mean_annual      = QS_mean_annual;
OUT.QF_mean_annual      = QF_mean_annual;
OUT.QT_mean_annual      = QT_mean_annual;
OUT.IE_mean_annual      = IE_mean_annual;
OUT.SE_mean_annual      = SE_mean_annual;
OUT.R_mean_annual       = R_mean_annual;
OUT.P_mean_annual       = P_mean_annual;
OUT.PET_mean_annual     = PET_mean_annual;
OUT.Ei_mean_annual      = Ei_mean_annual;
OUT.Et_mean_annual      = Et_mean_annual;
OUT.S_canopy_mean_annual = S_canopy_mean_annual;
OUT.pot_inf_mean_annual  = pot_inf_mean_annual;

end