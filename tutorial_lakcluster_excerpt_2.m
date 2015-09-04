clear ZRa_ ZCa_ ZRz_ ZCz_ ;
ZRa_ = cell(ncovs);ZCa_ = cell(ncovs);
ZRz_ = cell(ncovs);ZCz_ = cell(ncovs);
cov_bother_flag_AAAA=[]; cov_bother_flag_AZZA=[];
for ncov=0:ncovs-1;
lt__on = length(cov_A_up_{1+ncov}); lt_off = length(cov_A_dn_{1+ncov}); 
ls__on = length(cov_Z_up_{1+ncov}); ls_off = length(cov_Z_dn_{1+ncov});
cov_bother_flag_AAAA(1+ncov) = (lt__on>2) & (lt_off>2) ;
cov_bother_flag_AZZA(1+ncov) = (ls__on>2) & (ls_off>2) ;
if (cov_bother_flag_AAAA(1+ncov));
MA = lt__on + lt_off; NA = ncols_rem;
[tmp,cov_A_up_ij,tmp] = intersect(transpose(1:MA),cov_A_up_{1+ncov},'stable');
[tmp,cov_A_dn_ij,tmp] = intersect(transpose(1:MA),cov_A_dn_{1+ncov},'stable');
ZRa__on__on = XR_AAAA_upup_{1+ncov} - NA*(NA + lt__on - 1);
ZCa__on__on = XC_AAAA_upup_{1+ncov} - lt__on*(lt__on + NA - 1);
ZRa__on_off = XR_AAAA_updn_{1+ncov} - NA*lt_off;
ZCa__on_off = XC_AAAA_updn_{1+ncov} - lt__on*lt_off;
ZRa_off__on = XR_AAAA_dnup_{1+ncov} - NA*lt__on;
ZCa_off__on = XC_AAAA_dnup_{1+ncov} - lt_off*lt__on;
ZRa_off_off = XR_AAAA_dndn_{1+ncov} - NA*(NA + lt_off - 1);
ZCa_off_off = XC_AAAA_dndn_{1+ncov} - lt_off*(lt_off + NA - 1);
pa = (lt__on)/(lt__on+lt_off); qa = (lt_off)/(lt__on+lt_off);
ZRa__on = ...
        + min(ZRa__on__on*pa/pa,ZRa__on_off*pa/qa) ...
        + min(ZRa__on__on*qa/pa,ZRa__on_off*qa/qa) ...
  ;
ZRa_off = ...
        + min(ZRa_off__on*pa/pa,ZRa_off_off*pa/qa) ... 
        + min(ZRa_off__on*qa/pa,ZRa_off_off*qa/qa) ...
  ;
ZRa_{1+ncov}(cov_A_up_ij) = ZRa__on;
ZRa_{1+ncov}(cov_A_dn_ij) = ZRa_off;
papa = lt__on/(lt__on+lt_off)*lt__on/(lt__on+lt_off); 
paqa = lt__on/(lt__on+lt_off)*lt_off/(lt__on+lt_off);
qapa = lt_off/(lt__on+lt_off)*lt__on/(lt__on+lt_off); 
qaqa = lt_off/(lt__on+lt_off)*lt_off/(lt__on+lt_off);
ZCa_{1+ncov} = ...
  + min( ZCa__on__on*papa/papa , min( ZCa__on_off*papa/paqa , min( ZCa_off__on*papa/qapa , ZCa_off_off*papa/qaqa ))) ...
  + min( ZCa__on__on*paqa/papa , min( ZCa__on_off*paqa/paqa , min( ZCa_off__on*paqa/qapa , ZCa_off_off*paqa/qaqa ))) ...
  + min( ZCa__on__on*qapa/papa , min( ZCa__on_off*qapa/paqa , min( ZCa_off__on*qapa/qapa , ZCa_off_off*qapa/qaqa ))) ...
  + min( ZCa__on__on*qaqa/papa , min( ZCa__on_off*qaqa/paqa , min( ZCa_off__on*qaqa/qapa , ZCa_off_off*qaqa/qaqa ))) ...
  ;
ZRa_{1+ncov} = ZRa_{1+ncov} * (MA/MA)*(NA/NA)*(NA/NA);
ZCa_{1+ncov} = ZCa_{1+ncov} * (MA/MA)*(MA/MA)*(NA/NA);
end;%if (cov_bother_flag_AAAA(1+ncov));
MA = lt__on + lt_off; MZ = ls__on + ls_off; NA = ncols_rem; NZ = ncols_rem;
[tmp,cov_A_up_ij,tmp] = intersect(transpose(1:MA),cov_A_up_{1+ncov},'stable');
[tmp,cov_A_dn_ij,tmp] = intersect(transpose(1:MA),cov_A_dn_{1+ncov},'stable');
ZRz__on__on = XR_AZZA_upup_{1+ncov} - NA*ls__on;
ZCz__on__on = XC_AZZA_upup_{1+ncov} - lt__on*ls__on;
ZRz__on_off = XR_AZZA_updn_{1+ncov} - NA*ls_off;
ZCz__on_off = XC_AZZA_updn_{1+ncov} - lt__on*ls_off;
ZRz_off__on = XR_AZZA_dnup_{1+ncov} - NA*ls__on;
ZCz_off__on = XC_AZZA_dnup_{1+ncov} - lt_off*ls__on;
ZRz_off_off = XR_AZZA_dndn_{1+ncov} - NA*ls_off;
ZCz_off_off = XC_AZZA_dndn_{1+ncov} - lt_off*ls_off;
if ( cov_bother_flag_AAAA(1+ncov) &  cov_bother_flag_AZZA(1+ncov));
pz = (ls__on)/(ls__on+ls_off); qz = (ls_off)/(ls__on+ls_off);
ZRz__on = ...
        + min(ZRz__on__on*pz/pz,ZRz__on_off*pz/qz) ...
        + min(ZRz__on__on*qz/pz,ZRz__on_off*qz/qz) ...
  ;
ZRz_off = ...
        + min(ZRz_off__on*pz/pz,ZRz_off_off*pz/qz) ... 
        + min(ZRz_off__on*qz/pz,ZRz_off_off*qz/qz) ...
  ;
ZRz_{1+ncov}(cov_A_up_ij) = ZRz__on;
ZRz_{1+ncov}(cov_A_dn_ij) = ZRz_off;
papz = lt__on/(lt__on+lt_off)*ls__on/(ls__on+ls_off); 
paqz = lt__on/(lt__on+lt_off)*ls_off/(ls__on+ls_off);
qapz = lt_off/(lt__on+lt_off)*ls__on/(ls__on+ls_off); 
qaqz = lt_off/(lt__on+lt_off)*ls_off/(ls__on+ls_off);
ZCz_{1+ncov} = ...
  + min( ZCz__on__on*papz/papz , min( ZCz__on_off*papz/paqz , min( ZCz_off__on*papz/qapz , ZCz_off_off*papz/qaqz ))) ...
  + min( ZCz__on__on*paqz/papz , min( ZCz__on_off*paqz/paqz , min( ZCz_off__on*paqz/qapz , ZCz_off_off*paqz/qaqz ))) ...
  + min( ZCz__on__on*qapz/papz , min( ZCz__on_off*qapz/paqz , min( ZCz_off__on*qapz/qapz , ZCz_off_off*qapz/qaqz ))) ...
  + min( ZCz__on__on*qaqz/papz , min( ZCz__on_off*qaqz/paqz , min( ZCz_off__on*qaqz/qapz , ZCz_off_off*qaqz/qaqz ))) ...
  ;
ZRz_{1+ncov} = ZRz_{1+ncov} * (MA/MZ)*(NA/NZ)*(NA/NZ);
ZCz_{1+ncov} = ZCz_{1+ncov} * (MA/MZ)*(MA/MZ)*(NA/NZ);
end;%if ( cov_bother_flag_AAAA(1+ncov) &  cov_bother_flag_AZZA(1+ncov));
if (~cov_bother_flag_AAAA(1+ncov) &  cov_bother_flag_AZZA(1+ncov));
pz = (ls__on)/(ls__on+ls_off); qz = (ls_off)/(ls__on+ls_off);
ZRz__on = ...
        + min(ZRz__on__on*pz/pz,ZRz__on_off*pz/qz) ...
        + min(ZRz__on__on*qz/pz,ZRz__on_off*qz/qz) ...
  ;
ZRz_off = ...
        + min(ZRz_off__on*pz/pz,ZRz_off_off*pz/qz) ... 
        + min(ZRz_off__on*qz/pz,ZRz_off_off*qz/qz) ...
  ;
ZRz_{1+ncov}(cov_A_up_ij) = ZRz__on;
ZRz_{1+ncov}(cov_A_dn_ij) = ZRz_off;
papz = lt__on/(lt__on+lt_off)*ls__on/(ls__on+ls_off); 
paqz = lt__on/(lt__on+lt_off)*ls_off/(ls__on+ls_off);
qapz = lt_off/(lt__on+lt_off)*ls__on/(ls__on+ls_off); 
qaqz = lt_off/(lt__on+lt_off)*ls_off/(ls__on+ls_off);
if (lt__on>0);
ZCz_{1+ncov} = ...
  + min( ZCz__on__on*papz/papz , ZCz__on_off*papz/paqz ) ...
  + min( ZCz__on__on*paqz/papz , ZCz__on_off*paqz/paqz ) ...
  + min( ZCz__on__on*qapz/papz , ZCz__on_off*qapz/paqz ) ...
  + min( ZCz__on__on*qaqz/papz , ZCz__on_off*qaqz/paqz ) ...
  ;
end;%if (lt__on>0);
if (lt_off>0);
ZCz_{1+ncov} = ...
  + min( ZCz_off__on*papz/qapz , ZCz_off_off*papz/qaqz ) ...
  + min( ZCz_off__on*paqz/qapz , ZCz_off_off*paqz/qaqz ) ...
  + min( ZCz_off__on*qapz/qapz , ZCz_off_off*qapz/qaqz ) ...
  + min( ZCz_off__on*qaqz/qapz , ZCz_off_off*qaqz/qaqz ) ...
  ;
end;%if (lt_off>0);
ZRz_{1+ncov} = ZRz_{1+ncov} * (MA/MZ)*(NA/NZ)*(NA/NZ);
ZCz_{1+ncov} = ZCz_{1+ncov} * (MA/MZ)*(MA/MZ)*(NA/NZ);
end;%if (~cov_bother_flag_AAAA(1+ncov) &  cov_bother_flag_AZZA(1+ncov));
if ( cov_bother_flag_AAAA(1+ncov) & ~cov_bother_flag_AZZA(1+ncov));
pz = (ls__on)/(ls__on+ls_off); qz = (ls_off)/(ls__on+ls_off);
if (ls__on>0);
ZRz__on = ...
        + ZRz__on__on*pz/pz ...
        + ZRz__on__on*qz/pz ...
  ;
ZRz_off = ...
        + ZRz_off__on*pz/pz ... 
        + ZRz_off__on*qz/pz ...
  ;
end;%if (ls__on>0);
if (ls_off>0);
ZRz__on = ...
        + ZRz__on_off*pz/qz ...
        + ZRz__on_off*qz/qz ...
  ;
ZRz_off = ...
        + ZRz_off_off*pz/qz ... 
        + ZRz_off_off*qz/qz ...
  ;
end;%if (ls_off>0);
ZRz_{1+ncov}(cov_A_up_ij) = ZRz__on;
ZRz_{1+ncov}(cov_A_dn_ij) = ZRz_off;
papz = lt__on/(lt__on+lt_off)*ls__on/(ls__on+ls_off); 
paqz = lt__on/(lt__on+lt_off)*ls_off/(ls__on+ls_off);
qapz = lt_off/(lt__on+lt_off)*ls__on/(ls__on+ls_off); 
qaqz = lt_off/(lt__on+lt_off)*ls_off/(ls__on+ls_off);
if (ls__on>0);
ZCz_{1+ncov} = ...
  + min( ZCz__on__on*papz/papz , ZCz_off__on*papz/qapz ) ...
  + min( ZCz__on__on*paqz/papz , ZCz_off__on*paqz/qapz ) ...
  + min( ZCz__on__on*qapz/papz , ZCz_off__on*qapz/qapz ) ...
  + min( ZCz__on__on*qaqz/papz , ZCz_off__on*qaqz/qapz ) ...
  ;
end;%if (ls_on>0);
if (ls_off>0);
ZCz_{1+ncov} = ...
  + min( ZCz__on_off*papz/paqz , ZCz_off_off*papz/qaqz ) ...
  + min( ZCz__on_off*paqz/paqz , ZCz_off_off*paqz/qaqz ) ...
  + min( ZCz__on_off*qapz/paqz , ZCz_off_off*qapz/qaqz ) ...
  + min( ZCz__on_off*qaqz/paqz , ZCz_off_off*qaqz/qaqz ) ...
  ;
end;%if (ls_off>0);
ZRz_{1+ncov} = ZRz_{1+ncov} * (MA/MZ)*(NA/NZ)*(NA/NZ);
ZCz_{1+ncov} = ZCz_{1+ncov} * (MA/MZ)*(MA/MZ)*(NA/NZ);
end;%if (~cov_bother_flag_AAAA(1+ncov) & cov_bother_flag_AZZA(1+ncov));
end;%for ncov=0:ncovs-1;
cov_set_AAAA = 0;
for ncov=0:ncovs-1;
if (cov_bother_flag_AAAA(1+ncov) & ~cov_set_AAAA);
ZRa_min = ZRa_{1+ncov} ; ZCa_min = ZCa_{1+ncov} ; cov_set_AAAA=1;
elseif (cov_bother_flag_AAAA(1+ncov) & cov_set_AAAA);
ZRa_min = min(ZRa_min,ZRa_{1+ncov}) ; ZCa_min = min(ZCa_min,ZCa_{1+ncov}) ; 
end;% if cov_set_AAAA;
end;%for ncov=0:ncovs-1;
if ~cov_set_AAAA;
MA = nrows_rem; NA = ncols_rem;
ZRa_min = XR_AAAA_full - NA*(NA + MA - 1);
ZCa_min = XC_AAAA_full - MA*(MA + NA - 1);
ZRa_min = ZRa_min * (MA/MA)*(NA/NA)*(NA/NA);
ZCa_min = ZCa_min * (MA/MA)*(MA/MA)*(NA/NA);
cov_set_AAAA=1;
end;%if ~cov_set_AAAA;
cov_set_AZZA = 0;
for ncov=0:ncovs-1;
if ((cov_bother_flag_AAAA(1+ncov) | cov_bother_flag_AZZA(1+ncov)) & ~cov_set_AZZA);
ZRz_min = ZRz_{1+ncov}; ZCz_min = ZCz_{1+ncov} ; cov_set_AZZA=1;
elseif ((cov_bother_flag_AAAA(1+ncov) | cov_bother_flag_AZZA(1+ncov)) & cov_set_AZZA);
ZRz_min = min(ZRz_min,ZRz_{1+ncov}); ZCz_min = min(ZCz_min,ZCz_{1+ncov}) ; 
end;% if cov_set_AAAA | cov_set_AZZA;
end;%for ncov=0:ncovs-1;
if ~cov_set_AZZA;
MA = nrows_rem; NA = ncols_rem;
MZ = size(Z_orig,1); NZ = ncols_rem;
ZRz_min = XR_AZZA_full - NA*MZ;
ZCz_min = XC_AZZA_full - MA*MZ;
ZRz_min = ZRz_min * (MA/MZ)*(NA/NZ)*(NA/NZ);
ZCz_min = ZCz_min * (MA/MZ)*(MA/MZ)*(NA/NZ);
cov_set_AZZA=1;
end;%if ~cov_set_AZZA;
ZR = ZRa_min(:) - ZRz_min(:) ; ZC = ZCa_min(:) - ZCz_min(:) ;
if test_flag; disp(sprintf(' %% testing: ZR error = exp(%0.2f); ZC error = exp(%0.2f)',log(norm(ZR-ZR_test)),log(norm(ZC-ZC_test)))); end;%if test_flag;
ncovs_rem = sum(cov_bother_flag_AAAA)+sum(cov_bother_flag_AZZA);
