clear ZRz_ ZRz_g_ ZRz_l_ ZCz_ ZCz_g_ ZCz_l_ ;
ZRz_ = cell(ncovs);ZCz_ = cell(ncovs);
ZRz_g_ = cell(ncovs);ZCz_g_ = cell(ncovs);
ZRz_l_ = cell(ncovs);ZCz_l_ = cell(ncovs);
cov_bother_flag_AAAA=[]; cov_bother_flag_AZZA=[];
for ncov=1:ncovs;
lt__on = length(cov_A_up_{ncov}); lt_off = length(cov_A_dn_{ncov}); 
ls__on = length(cov_Z_up_{ncov}); ls_off = length(cov_Z_dn_{ncov});
cov_bother_flag_AAAA(ncov) = (lt__on>2) & (lt_off>2) ;
cov_bother_flag_AZZA(ncov) = (ls__on>2) & (ls_off>2) ;
MA = lt__on + lt_off; MZ = ls__on + ls_off; NA = ncols_rem; NZ = ncols_rem;
ZRz_g_xx_up = X_P_xx_up{ncov}; ZRz_l_xx_up = ls__on - X_P_xx_up{ncov};
ZRz_g_xx_dn = X_P_xx_dn{ncov}; ZRz_l_xx_dn = ls_off - X_P_xx_dn{ncov};
ZCz_g_up_up = sum(X_P_up_up{ncov},1); ZCz_l_up_up = lt__on*ls__on - ZCz_g_up_up;
ZCz_g_up_dn = sum(X_P_up_dn{ncov},1); ZCz_l_up_dn = lt__on*ls_off - ZCz_g_up_dn;
ZCz_g_dn_up = sum(X_P_dn_up{ncov},1); ZCz_l_dn_up = lt_off*ls__on - ZCz_g_dn_up;
ZCz_g_dn_dn = sum(X_P_dn_dn{ncov},1); ZCz_l_dn_dn = lt_off*ls_off - ZCz_g_dn_dn;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( cov_bother_flag_AAAA(ncov) &  cov_bother_flag_AZZA(ncov));
pz = (ls__on)/(ls__on+ls_off); qz = (ls_off)/(ls__on+ls_off);
ZRz_g_{ncov} = ...
        + min(ZRz_g_xx_up*pz/pz,ZRz_g_xx_dn*pz/qz) ...
        + min(ZRz_g_xx_up*qz/pz,ZRz_g_xx_dn*qz/qz) ...
  ;
ZRz_l_{ncov} = ...
        + min(ZRz_l_xx_up*pz/pz,ZRz_l_xx_dn*pz/qz) ...
        + min(ZRz_l_xx_up*qz/pz,ZRz_l_xx_dn*qz/qz) ...
  ;
ZRz_{ncov}(r_rem{iteration},c_rem{iteration}) = max(ZRz_g_{ncov}(r_rem{iteration},c_rem{iteration}),ZRz_l_{ncov}(r_rem{iteration},c_rem{iteration})) ;
papz = lt__on/(lt__on+lt_off)*ls__on/(ls__on+ls_off); 
paqz = lt__on/(lt__on+lt_off)*ls_off/(ls__on+ls_off);
qapz = lt_off/(lt__on+lt_off)*ls__on/(ls__on+ls_off); 
qaqz = lt_off/(lt__on+lt_off)*ls_off/(ls__on+ls_off);
ZCz_g_{ncov} = ...
  + min( ZCz_g_up_up*papz/papz , min( ZCz_g_up_dn*papz/paqz , min( ZCz_g_dn_up*papz/qapz , ZCz_g_dn_dn*papz/qaqz ))) ...
  + min( ZCz_g_up_up*paqz/papz , min( ZCz_g_up_dn*paqz/paqz , min( ZCz_g_dn_up*paqz/qapz , ZCz_g_dn_dn*paqz/qaqz ))) ...
  + min( ZCz_g_up_up*qapz/papz , min( ZCz_g_up_dn*qapz/paqz , min( ZCz_g_dn_up*qapz/qapz , ZCz_g_dn_dn*qapz/qaqz ))) ...
  + min( ZCz_g_up_up*qaqz/papz , min( ZCz_g_up_dn*qaqz/paqz , min( ZCz_g_dn_up*qaqz/qapz , ZCz_g_dn_dn*qaqz/qaqz ))) ...
  ;
ZCz_l_{ncov} = ...
  + min( ZCz_l_up_up*papz/papz , min( ZCz_l_up_dn*papz/paqz , min( ZCz_l_dn_up*papz/qapz , ZCz_l_dn_dn*papz/qaqz ))) ...
  + min( ZCz_l_up_up*paqz/papz , min( ZCz_l_up_dn*paqz/paqz , min( ZCz_l_dn_up*paqz/qapz , ZCz_l_dn_dn*paqz/qaqz ))) ...
  + min( ZCz_l_up_up*qapz/papz , min( ZCz_l_up_dn*qapz/paqz , min( ZCz_l_dn_up*qapz/qapz , ZCz_l_dn_dn*qapz/qaqz ))) ...
  + min( ZCz_l_up_up*qaqz/papz , min( ZCz_l_up_dn*qaqz/paqz , min( ZCz_l_dn_up*qaqz/qapz , ZCz_l_dn_dn*qaqz/qaqz ))) ...
  ;
ZCz_{ncov}(1,c_rem{iteration}) = max(ZCz_g_{ncov}(1,c_rem{iteration}),ZCz_l_{ncov}(1,c_rem{iteration})) ;
end;%if ( cov_bother_flag_AAAA(ncov) &  cov_bother_flag_AZZA(ncov));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (~cov_bother_flag_AAAA(ncov) &  cov_bother_flag_AZZA(ncov));
pz = (ls__on)/(ls__on+ls_off); qz = (ls_off)/(ls__on+ls_off);
ZRz_g_{ncov} = ...
        + min(ZRz_g_xx_up*pz/pz,ZRz_g_xx_dn*pz/qz) ...
        + min(ZRz_g_xx_up*qz/pz,ZRz_g_xx_dn*qz/qz) ...
  ;
ZRz_l_{ncov} = ...
        + min(ZRz_l_xx_up*pz/pz,ZRz_l_xx_dn*pz/qz) ...
        + min(ZRz_l_xx_up*qz/pz,ZRz_l_xx_dn*qz/qz) ...
  ;
ZRz_{ncov}(r_rem{iteration},c_rem{iteration}) = max(ZRz_g_{ncov}(r_rem{iteration},c_rem{iteration}),ZRz_l_{ncov}(r_rem{iteration},c_rem{iteration})) ;
papz = lt__on/(lt__on+lt_off)*ls__on/(ls__on+ls_off); 
paqz = lt__on/(lt__on+lt_off)*ls_off/(ls__on+ls_off);
qapz = lt_off/(lt__on+lt_off)*ls__on/(ls__on+ls_off); 
qaqz = lt_off/(lt__on+lt_off)*ls_off/(ls__on+ls_off);
if (lt__on>0);
ZCz_g_{ncov} = ...
  + min( ZCz_g_up_up*papz/papz , ZCz_g_up_dn*papz/paqz ) ...
  + min( ZCz_g_up_up*paqz/papz , ZCz_g_up_dn*paqz/paqz ) ...
  + min( ZCz_g_up_up*qapz/papz , ZCz_g_up_dn*qapz/paqz ) ...
  + min( ZCz_g_up_up*qaqz/papz , ZCz_g_up_dn*qaqz/paqz ) ...
  ;
ZCz_l_{ncov} = ...
  + min( ZCz_l_up_up*papz/papz , ZCz_l_up_dn*papz/paqz ) ...
  + min( ZCz_l_up_up*paqz/papz , ZCz_l_up_dn*paqz/paqz ) ...
  + min( ZCz_l_up_up*qapz/papz , ZCz_l_up_dn*qapz/paqz ) ...
  + min( ZCz_l_up_up*qaqz/papz , ZCz_l_up_dn*qaqz/paqz ) ...
  ;
ZCz_{ncov}(1,c_rem{iteration}) = max(ZCz_g_{ncov}(1,c_rem{iteration}),ZCz_l_{ncov}(1,c_rem{iteration})) ;
end;%if (lt__on>0);
if (lt_off>0);
ZCz_g_{ncov} = ...
  + min( ZCz_g_dn_up*papz/qapz , ZCz_g_dn_dn*papz/qaqz ) ...
  + min( ZCz_g_dn_up*paqz/qapz , ZCz_g_dn_dn*paqz/qaqz ) ...
  + min( ZCz_g_dn_up*qapz/qapz , ZCz_g_dn_dn*qapz/qaqz ) ...
  + min( ZCz_g_dn_up*qaqz/qapz , ZCz_g_dn_dn*qaqz/qaqz ) ...
  ;
ZCz_l_{ncov} = ...
  + min( ZCz_l_dn_up*papz/qapz , ZCz_l_dn_dn*papz/qaqz ) ...
  + min( ZCz_l_dn_up*paqz/qapz , ZCz_l_dn_dn*paqz/qaqz ) ...
  + min( ZCz_l_dn_up*qapz/qapz , ZCz_l_dn_dn*qapz/qaqz ) ...
  + min( ZCz_l_dn_up*qaqz/qapz , ZCz_l_dn_dn*qaqz/qaqz ) ...
  ;
ZCz_{ncov}(1,c_rem{iteration}) = max(ZCz_g_{ncov}(1,c_rem{iteration}),ZCz_l_{ncov}(1,c_rem{iteration})) ;
end;%if (lt_off>0);
end;%if (~cov_bother_flag_AAAA(ncov) &  cov_bother_flag_AZZA(ncov));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( cov_bother_flag_AAAA(ncov) & ~cov_bother_flag_AZZA(ncov));
pz = (ls__on)/(ls__on+ls_off); qz = (ls_off)/(ls__on+ls_off);
if (ls__on>0);
ZRz_g_{ncov} = ...
        + (ZRz_g_xx_up*pz/pz) ...
        + (ZRz_g_xx_up*qz/pz) ...
  ;
ZRz_l_{ncov} = ...
        + (ZRz_l_xx_up*pz/pz) ...
        + (ZRz_l_xx_up*qz/pz) ...
  ;
ZRz_{ncov}(r_rem{iteration},c_rem{iteration}) = max(ZRz_g_{ncov}(r_rem{iteration},c_rem{iteration}),ZRz_l_{ncov}(r_rem{iteration},c_rem{iteration})) ;
end;%if (ls__on>0);
if (ls_off>0);
ZRz_g_{ncov} = ...
        + (ZRz_g_xx_dn*pz/qz) ...
        + (ZRz_g_xx_dn*qz/qz) ...
  ;
ZRz_l_{ncov} = ...
        + (ZRz_l_xx_dn*pz/qz) ...
        + (ZRz_l_xx_dn*qz/qz) ...
  ;
ZRz_{ncov}(r_rem{iteration},c_rem{iteration}) = max(ZRz_g_{ncov}(r_rem{iteration},c_rem{iteration}),ZRz_l_{ncov}(r_rem{iteration},c_rem{iteration})) ;
end;%if (ls_off>0);
papz = lt__on/(lt__on+lt_off)*ls__on/(ls__on+ls_off); 
paqz = lt__on/(lt__on+lt_off)*ls_off/(ls__on+ls_off);
qapz = lt_off/(lt__on+lt_off)*ls__on/(ls__on+ls_off); 
qaqz = lt_off/(lt__on+lt_off)*ls_off/(ls__on+ls_off);
if (ls__on>0);
ZCz_g_{ncov} = ...
  + min( ZCz_g_up_up*papz/papz , ZCz_g_dn_up*papz/qapz ) ...
  + min( ZCz_g_up_up*paqz/papz , ZCz_g_dn_up*paqz/qapz ) ...
  + min( ZCz_g_up_up*qapz/papz , ZCz_g_dn_up*qapz/qapz ) ...
  + min( ZCz_g_up_up*qaqz/papz , ZCz_g_dn_up*qaqz/qapz ) ...
  ;
ZCz_l_{ncov} = ...
  + min( ZCz_l_up_up*papz/papz , ZCz_l_dn_up*papz/qapz ) ...
  + min( ZCz_l_up_up*paqz/papz , ZCz_l_dn_up*paqz/qapz ) ...
  + min( ZCz_l_up_up*qapz/papz , ZCz_l_dn_up*qapz/qapz ) ...
  + min( ZCz_l_up_up*qaqz/papz , ZCz_l_dn_up*qaqz/qapz ) ...
  ;
ZCz_{ncov}(1,c_rem{iteration}) = max(ZCz_g_{ncov}(1,c_rem{iteration}),ZCz_l_{ncov}(1,c_rem{iteration})) ;
end;%if (ls_on>0);
if (ls_off>0);
ZCz_g_{ncov} = ...
  + min( ZCz_g_up_dn*papz/paqz , ZCz_g_dn_dn*papz/qaqz ) ...
  + min( ZCz_g_up_dn*paqz/paqz , ZCz_g_dn_dn*paqz/qaqz ) ...
  + min( ZCz_g_up_dn*qapz/paqz , ZCz_g_dn_dn*qapz/qaqz ) ...
  + min( ZCz_g_up_dn*qaqz/paqz , ZCz_g_dn_dn*qaqz/qaqz ) ...
  ;
ZCz_l_{ncov} = ...
  + min( ZCz_l_up_dn*papz/paqz , ZCz_l_dn_dn*papz/qaqz ) ...
  + min( ZCz_l_up_dn*paqz/paqz , ZCz_l_dn_dn*paqz/qaqz ) ...
  + min( ZCz_l_up_dn*qapz/paqz , ZCz_l_dn_dn*qapz/qaqz ) ...
  + min( ZCz_l_up_dn*qaqz/paqz , ZCz_l_dn_dn*qaqz/qaqz ) ...
  ;
ZCz_{ncov}(1,c_rem{iteration}) = max(ZCz_g_{ncov}(1,c_rem{iteration}),ZCz_l_{ncov}(1,c_rem{iteration})) ;
end;%if (ls_off>0);
end;%if (~cov_bother_flag_AAAA(ncov) & cov_bother_flag_AZZA(ncov));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end;%for ncov=1:ncovs;
cov_set_AZZA = 0;
for ncov=1:ncovs;
if ((cov_bother_flag_AAAA(ncov) | cov_bother_flag_AZZA(ncov)) & ~cov_set_AZZA);
ZRz_min = ZRz_{ncov}; ZCz_min = ZCz_{ncov} ; cov_set_AZZA=1;
elseif ((cov_bother_flag_AAAA(ncov) | cov_bother_flag_AZZA(ncov)) & cov_set_AZZA);
ZRz_min = min(ZRz_min,ZRz_{ncov}); ZCz_min = min(ZCz_min,ZCz_{ncov}) ; 
end;% if cov_set_AAAA | cov_set_AZZA;
end;%for ncov=1:ncovs;
if ~cov_set_AZZA;
MA = nrows_rem; NA = ncols_rem; MZ = size(Z_orig,1); NZ = ncols_rem;
ZRz_min = max(X_P_all,MZ-X_P_all);
tmp = sum(X_P_all,1);
ZCz_min = max(tmp,MA*MZ-tmp);
cov_set_AZZA=1;
end;%if ~cov_set_AZZA;
ZR = sum(ZRz_min(1:nrows_rem,1:ncols_rem),2); ZR = ZR(:); 
ZC = ZCz_min(:); ZC = ZC(1:ncols_rem);
