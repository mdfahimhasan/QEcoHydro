function [Cl_Q,Cl_ET] = get_conc_out(Cl_S,OMEGA_Q_SAVE,OMEGA_ET_SAVE)

clear mass_cons Cl_Q Cl_ET

for i = 3000:-1:2
             
    omega_Q       = [OMEGA_Q_SAVE(i,1);  diff( OMEGA_Q_SAVE(:,i) )  ] ;
    omega_ET      = [OMEGA_ET_SAVE(i,1);  diff( OMEGA_ET_SAVE(:,i) )  ] ;

    Cl_history   = Cl_S(:,i-1);

    %%%% Make sure mass is being converved
    vec_aux           = Cl_history./Cl_history;
    test              = omega_Q.*vec_aux;
    mass_cons(i,1)    = nansum(test)-nansum(omega_Q);
    %%%% 

    
    Cl_Q(i,1)          = nansum( omega_Q./nansum(omega_Q).*Cl_history );
    Cl_ET(i,1)         = nansum( omega_ET./nansum(omega_ET).*Cl_history );

end