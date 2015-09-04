if (test_flag);
clear ZRa_ ZCa_ ZRz_ ZCz_ ;
ZRa_ = cell(ncovs);ZCa_ = cell(ncovs);
ZRz_ = cell(ncovs);ZCz_ = cell(ncovs);
for ncov=0:ncovs-1;
T__on_{1+ncov} = find(T(:,1+1+ncov)==1);T_off_{1+ncov} = find(T(:,1+1+ncov)==0);
S__on_{1+ncov} = find(S(:,1+1+ncov)==1);S_off_{1+ncov} = find(S(:,1+1+ncov)==0);
lt__on = length(T__on_{1+ncov}); lt_off = length(T_off_{1+ncov}); 
ls__on = length(S__on_{1+ncov}); ls_off = length(S_off_{1+ncov});
cov_bother_flag_AAAA(1+ncov) = (lt__on>2) & (lt_off>2) ;
cov_bother_flag_AZZA(1+ncov) = (ls__on>2) & (ls_off>2) ;
if (cov_bother_flag_AAAA(1+ncov));
MA = lt__on + lt_off; NA = size(A,2);
[tmp,T__on_ij,tmp] = intersect(transpose(1:MA),T__on_{1+ncov},'stable');
[tmp,T_off_ij,tmp] = intersect(transpose(1:MA),T_off_{1+ncov},'stable');
tmp__on__on = (A(T__on_{1+ncov},:)*transpose(A(T__on_{1+ncov},:)));
ZRa__on__on = diag((tmp__on__on)*transpose(tmp__on__on)) - NA*(NA + lt__on - 1);
ZCa__on__on = diag(transpose(A(T__on_{1+ncov},:))*(tmp__on__on)*A(T__on_{1+ncov},:)) - lt__on*(lt__on + NA - 1);
tmp__on_off = (A(T__on_{1+ncov},:)*transpose(A(T_off_{1+ncov},:)));
ZRa__on_off = diag((tmp__on_off)*transpose(tmp__on_off)) - NA*lt_off;
ZCa__on_off = diag(transpose(A(T__on_{1+ncov},:))*(tmp__on_off)*A(T_off_{1+ncov},:)) - lt__on*lt_off;
tmp_off__on = (A(T_off_{1+ncov},:)*transpose(A(T__on_{1+ncov},:)));
ZRa_off__on = diag((tmp_off__on)*transpose(tmp_off__on)) - NA*lt__on;
ZCa_off__on = diag(transpose(A(T_off_{1+ncov},:))*(tmp_off__on)*A(T__on_{1+ncov},:)) - lt_off*lt__on;
tmp_off_off = (A(T_off_{1+ncov},:)*transpose(A(T_off_{1+ncov},:)));
ZRa_off_off = diag((tmp_off_off)*transpose(tmp_off_off)) - NA*(NA + lt_off - 1);
ZCa_off_off = diag(transpose(A(T_off_{1+ncov},:))*(tmp_off_off)*A(T_off_{1+ncov},:)) - lt_off*(lt_off + NA - 1);
pa = (lt__on)/(lt__on+lt_off); qa = (lt_off)/(lt__on+lt_off);
ZRa__on = ...
        + min(ZRa__on__on*pa/pa,ZRa__on_off*pa/qa) ...
        + min(ZRa__on__on*qa/pa,ZRa__on_off*qa/qa) ...
  ;
ZRa_off = ...
        + min(ZRa_off__on*pa/pa,ZRa_off_off*pa/qa) ... 
        + min(ZRa_off__on*qa/pa,ZRa_off_off*qa/qa) ...
  ;
ZRa_{1+ncov}(T__on_ij) = ZRa__on;
ZRa_{1+ncov}(T_off_ij) = ZRa_off;
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
if ( cov_bother_flag_AAAA(1+ncov) &  cov_bother_flag_AZZA(1+ncov));
MA = lt__on + lt_off; MZ = ls__on + ls_off; NA = size(A,2); NZ = size(Z,2);
[tmp,T__on_ij,tmp] = intersect(transpose(1:MA),T__on_{1+ncov},'stable');
[tmp,T_off_ij,tmp] = intersect(transpose(1:MA),T_off_{1+ncov},'stable');
tmp__on__on = (A(T__on_{1+ncov},:)*transpose(Z(S__on_{1+ncov},:)));
ZRz__on__on = diag((tmp__on__on)*transpose(tmp__on__on)) - NA*ls__on;
ZCz__on__on = diag(transpose(A(T__on_{1+ncov},:))*(tmp__on__on)*Z(S__on_{1+ncov},:)) - lt__on*ls__on;
tmp__on_off = (A(T__on_{1+ncov},:)*transpose(Z(S_off_{1+ncov},:)));
ZRz__on_off = diag((tmp__on_off)*transpose(tmp__on_off)) - NA*ls_off;
ZCz__on_off = diag(transpose(A(T__on_{1+ncov},:))*(tmp__on_off)*Z(S_off_{1+ncov},:)) - lt__on*ls_off;
tmp_off__on = (A(T_off_{1+ncov},:)*transpose(Z(S__on_{1+ncov},:)));
ZRz_off__on = diag((tmp_off__on)*transpose(tmp_off__on)) - NA*ls__on;
ZCz_off__on = diag(transpose(A(T_off_{1+ncov},:))*(tmp_off__on)*Z(S__on_{1+ncov},:)) - lt_off*ls__on;
tmp_off_off = (A(T_off_{1+ncov},:)*transpose(Z(S_off_{1+ncov},:)));
ZRz_off_off = diag((tmp_off_off)*transpose(tmp_off_off)) - NA*ls_off;
ZCz_off_off = diag(transpose(A(T_off_{1+ncov},:))*(tmp_off_off)*Z(S_off_{1+ncov},:)) - lt_off*ls_off;
pz = (ls__on)/(ls__on+ls_off); qz = (ls_off)/(ls__on+ls_off);
ZRz__on = ...
        + min(ZRz__on__on*pz/pz,ZRz__on_off*pz/qz) ...
        + min(ZRz__on__on*qz/pz,ZRz__on_off*qz/qz) ...
  ;
ZRz_off = ...
        + min(ZRz_off__on*pz/pz,ZRz_off_off*pz/qz) ... 
        + min(ZRz_off__on*qz/pz,ZRz_off_off*qz/qz) ...
  ;
ZRz_{1+ncov}(T__on_ij) = ZRz__on;
ZRz_{1+ncov}(T_off_ij) = ZRz_off;
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
MA = lt__on + lt_off; MZ = ls__on + ls_off; NA = size(A,2); NZ = size(Z,2);
[tmp,T__on_ij,tmp] = intersect(transpose(1:MA),T__on_{1+ncov},'stable');
[tmp,T_off_ij,tmp] = intersect(transpose(1:MA),T_off_{1+ncov},'stable');
tmp__on__on = (A(T__on_{1+ncov},:)*transpose(Z(S__on_{1+ncov},:)));
ZRz__on__on = diag((tmp__on__on)*transpose(tmp__on__on)) - NA*ls__on;
ZCz__on__on = diag(transpose(A(T__on_{1+ncov},:))*(tmp__on__on)*Z(S__on_{1+ncov},:)) - lt__on*ls__on;
tmp__on_off = (A(T__on_{1+ncov},:)*transpose(Z(S_off_{1+ncov},:)));
ZRz__on_off = diag((tmp__on_off)*transpose(tmp__on_off)) - NA*ls_off;
ZCz__on_off = diag(transpose(A(T__on_{1+ncov},:))*(tmp__on_off)*Z(S_off_{1+ncov},:)) - lt__on*ls_off;
tmp_off__on = (A(T_off_{1+ncov},:)*transpose(Z(S__on_{1+ncov},:)));
ZRz_off__on = diag((tmp_off__on)*transpose(tmp_off__on)) - NA*ls__on;
ZCz_off__on = diag(transpose(A(T_off_{1+ncov},:))*(tmp_off__on)*Z(S__on_{1+ncov},:)) - lt_off*ls__on;
tmp_off_off = (A(T_off_{1+ncov},:)*transpose(Z(S_off_{1+ncov},:)));
ZRz_off_off = diag((tmp_off_off)*transpose(tmp_off_off)) - NA*ls_off;
ZCz_off_off = diag(transpose(A(T_off_{1+ncov},:))*(tmp_off_off)*Z(S_off_{1+ncov},:)) - lt_off*ls_off;
pz = (ls__on)/(ls__on+ls_off); qz = (ls_off)/(ls__on+ls_off);
ZRz__on = ...
        + min(ZRz__on__on*pz/pz,ZRz__on_off*pz/qz) ...
        + min(ZRz__on__on*qz/pz,ZRz__on_off*qz/qz) ...
  ;
ZRz_off = ...
        + min(ZRz_off__on*pz/pz,ZRz_off_off*pz/qz) ... 
        + min(ZRz_off__on*qz/pz,ZRz_off_off*qz/qz) ...
  ;
ZRz_{1+ncov}(T__on_ij) = ZRz__on;
ZRz_{1+ncov}(T_off_ij) = ZRz_off;
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
MA = lt__on + lt_off; MZ = ls__on + ls_off; NA = size(A,2); NZ = size(Z,2);
[tmp,T__on_ij,tmp] = intersect(transpose(1:MA),T__on_{1+ncov},'stable');
[tmp,T_off_ij,tmp] = intersect(transpose(1:MA),T_off_{1+ncov},'stable');
tmp__on__on = (A(T__on_{1+ncov},:)*transpose(Z(S__on_{1+ncov},:)));
ZRz__on__on = diag((tmp__on__on)*transpose(tmp__on__on)) - NA*ls__on;
ZCz__on__on = diag(transpose(A(T__on_{1+ncov},:))*(tmp__on__on)*Z(S__on_{1+ncov},:)) - lt__on*ls__on;
tmp__on_off = (A(T__on_{1+ncov},:)*transpose(Z(S_off_{1+ncov},:)));
ZRz__on_off = diag((tmp__on_off)*transpose(tmp__on_off)) - NA*ls_off;
ZCz__on_off = diag(transpose(A(T__on_{1+ncov},:))*(tmp__on_off)*Z(S_off_{1+ncov},:)) - lt__on*ls_off;
tmp_off__on = (A(T_off_{1+ncov},:)*transpose(Z(S__on_{1+ncov},:)));
ZRz_off__on = diag((tmp_off__on)*transpose(tmp_off__on)) - NA*ls__on;
ZCz_off__on = diag(transpose(A(T_off_{1+ncov},:))*(tmp_off__on)*Z(S__on_{1+ncov},:)) - lt_off*ls__on;
tmp_off_off = (A(T_off_{1+ncov},:)*transpose(Z(S_off_{1+ncov},:)));
ZRz_off_off = diag((tmp_off_off)*transpose(tmp_off_off)) - NA*ls_off;
ZCz_off_off = diag(transpose(A(T_off_{1+ncov},:))*(tmp_off_off)*Z(S_off_{1+ncov},:)) - lt_off*ls_off;
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
ZRz_{1+ncov}(T__on_ij) = ZRz__on;
ZRz_{1+ncov}(T_off_ij) = ZRz_off;
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
end;% for ncov=0:ncovs-1;
cov_set_AAAA = 0;
for ncov=0:ncovs-1;
if (cov_bother_flag_AAAA(1+ncov) & ~cov_set_AAAA);
ZRa_min = ZRa_{1+ncov} ; ZCa_min = ZCa_{1+ncov} ; cov_set_AAAA=1;
elseif (cov_bother_flag_AAAA(1+ncov) & cov_set_AAAA);
ZRa_min = min(ZRa_min,ZRa_{1+ncov}) ; ZCa_min = min(ZCa_min,ZCa_{1+ncov}) ; 
end;% if cov_set_AAAA;
end;%for ncov=0:ncovs-1;
if ~cov_set_AAAA;
[MA,NA] = size(A);
tmp = (A*transpose(A));
ZRa_min = diag(tmp*transpose(tmp)) - NA*(NA + MA - 1);
ZCa_min = diag(transpose(A)*(tmp)*A) - MA*(MA + NA - 1);
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
[MA,NA] = size(A);
[MZ,NZ] = size(Z);
tmp = (A*transpose(Z));
ZRz_min = diag((tmp)*transpose(tmp)) - NA*MZ;
ZCz_min = diag(transpose(A)*(tmp)*Z) - MA*MZ;
ZRz_min = ZRz_min * (MA/MZ)*(NA/NZ)*(NA/NZ);
ZCz_min = ZCz_min * (MA/MZ)*(MA/MZ)*(NA/NZ);
cov_set_AZZA=1;
end;%if ~cov_set_AZZA;
ZR_test = ZRa_min(:) - ZRz_min(:) ; ZC_test = ZCa_min(:) - ZCz_min(:) ;
end;%if (test_flag);
