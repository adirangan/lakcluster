function [LR,LC,Lx,Lx_AAAA,Lx_AZZA] = cfilter_AZZA_Y_0_getscore(A_orig,Z_orig,T_orig,S_orig);
% calculates loop score (i.e., sum of score vector) ;
% returns LR: row score ;
% returns LC: col score ;
% returns Lx: number of covariates remaining ;

ncovs = size(T_orig,2)-1; % the first column of T,S should be all ones ;
%gamma = GLOBAL_CFILTER_AGGRESSIVE;
%disp(sprintf('beginning cfilter: gamma %f',gamma));
[MA,NA] = size(A_orig); [MZ,NZ] = size(Z_orig);
out_xdrop = zeros(MA+NA,2); out_xdrop_ij=1;
%out_trace = zeros(MA+NA,6); 
out_trace = zeros(1,6); 
row_ij = 1:MA; col_ij = 1:NA;
A = A_orig; Z = Z_orig; T = T_orig>0; S = S_orig>0;
rij = [];cij = [];
iteration=1;
r_rem{iteration} = row_ij;c_rem{iteration} = col_ij;
%while min(size(A))>2;
%disp(sprintf('iteration %d, A size %d-x-%d = exp(%0.2f)-x-exp(%0.2f) ratio %0.2f',iteration,size(A),log(size(A)),size(A,2)/size(A,1)));
ncovs = size(T,2)-1; % the first column of T,S should be all ones ;
cov_bother_flag_AAAA=zeros(1+ncovs,1);cov_bother_flag_AAAA(1)=1;
cov_bother_flag_AZZA=zeros(1+ncovs,1);cov_bother_flag_AZZA(1)=1;
clear ZRa ZCa ZRz ZCz ;
ZRa = cell(ncovs);ZCa = cell(ncovs);
ZRz = cell(ncovs);ZCz = cell(ncovs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through covariates 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ncov=0:ncovs-1;
T__on{1+ncov} = find(T(:,1+1+ncov)==1);T_off{1+ncov} = find(T(:,1+1+ncov)==0);
S__on{1+ncov} = find(S(:,1+1+ncov)==1);S_off{1+ncov} = find(S(:,1+1+ncov)==0);
lt__on = length(T__on{1+ncov}); lt_off = length(T_off{1+ncov}); 
ls__on = length(S__on{1+ncov}); ls_off = length(S_off{1+ncov});
cov_bother_flag_AAAA(1+ncov) = (lt__on>2) & (lt_off>2) ;
cov_bother_flag_AZZA(1+ncov) = (ls__on>2) & (ls_off>2) ;

if (cov_bother_flag_AAAA(1+ncov));
MA = lt__on + lt_off; NA = size(A,2);
[tmp,T__on_ij,tmp] = intersect(transpose(1:MA),T__on{1+ncov},'stable');
[tmp,T_off_ij,tmp] = intersect(transpose(1:MA),T_off{1+ncov},'stable');
tmp__on__on = (A(T__on{1+ncov},:)*transpose(A(T__on{1+ncov},:)));
ZRa__on__on = diag((tmp__on__on)*transpose(tmp__on__on)) - NA*(NA + lt__on - 1);
ZCa__on__on = diag(transpose(A(T__on{1+ncov},:))*(tmp__on__on)*A(T__on{1+ncov},:)) - lt__on*(lt__on + NA - 1);
tmp__on_off = (A(T__on{1+ncov},:)*transpose(A(T_off{1+ncov},:)));
ZRa__on_off = diag((tmp__on_off)*transpose(tmp__on_off)) - NA*lt_off;
ZCa__on_off = diag(transpose(A(T__on{1+ncov},:))*(tmp__on_off)*A(T_off{1+ncov},:)) - lt__on*lt_off;
tmp_off__on = (A(T_off{1+ncov},:)*transpose(A(T__on{1+ncov},:)));
ZRa_off__on = diag((tmp_off__on)*transpose(tmp_off__on)) - NA*lt__on;
ZCa_off__on = diag(transpose(A(T_off{1+ncov},:))*(tmp_off__on)*A(T__on{1+ncov},:)) - lt_off*lt__on;
tmp_off_off = (A(T_off{1+ncov},:)*transpose(A(T_off{1+ncov},:)));
ZRa_off_off = diag((tmp_off_off)*transpose(tmp_off_off)) - NA*(NA + lt_off - 1);
ZCa_off_off = diag(transpose(A(T_off{1+ncov},:))*(tmp_off_off)*A(T_off{1+ncov},:)) - lt_off*(lt_off + NA - 1);
pa = (lt__on)/(lt__on+lt_off); qa = (lt_off)/(lt__on+lt_off);
ZRa__on = ...
        + min(ZRa__on__on*pa/pa,ZRa__on_off*pa/qa) ...
        + min(ZRa__on__on*qa/pa,ZRa__on_off*qa/qa) ...
  ;
ZRa_off = ...
        + min(ZRa_off__on*pa/pa,ZRa_off_off*pa/qa) ... 
        + min(ZRa_off__on*qa/pa,ZRa_off_off*qa/qa) ...
  ;
ZRa{1+ncov}(T__on_ij) = ZRa__on;
ZRa{1+ncov}(T_off_ij) = ZRa_off;
papa = lt__on/(lt__on+lt_off)*lt__on/(lt__on+lt_off); 
paqa = lt__on/(lt__on+lt_off)*lt_off/(lt__on+lt_off);
qapa = lt_off/(lt__on+lt_off)*lt__on/(lt__on+lt_off); 
qaqa = lt_off/(lt__on+lt_off)*lt_off/(lt__on+lt_off);
ZCa{1+ncov} = ...
  + min( ZCa__on__on*papa/papa , min( ZCa__on_off*papa/paqa , min( ZCa_off__on*papa/qapa , ZCa_off_off*papa/qaqa ))) ...
  + min( ZCa__on__on*paqa/papa , min( ZCa__on_off*paqa/paqa , min( ZCa_off__on*paqa/qapa , ZCa_off_off*paqa/qaqa ))) ...
  + min( ZCa__on__on*qapa/papa , min( ZCa__on_off*qapa/paqa , min( ZCa_off__on*qapa/qapa , ZCa_off_off*qapa/qaqa ))) ...
  + min( ZCa__on__on*qaqa/papa , min( ZCa__on_off*qaqa/paqa , min( ZCa_off__on*qaqa/qapa , ZCa_off_off*qaqa/qaqa ))) ...
  ;
ZRa{1+ncov} = ZRa{1+ncov} * (MA/MA)*(NA/NA)*(NA/NA);
ZCa{1+ncov} = ZCa{1+ncov} * (MA/MA)*(MA/MA)*(NA/NA);
end;%if (cov_bother_flag_AAAA(1+ncov));

if ( cov_bother_flag_AAAA(1+ncov) &  cov_bother_flag_AZZA(1+ncov));
MA = lt__on + lt_off; MZ = ls__on + ls_off; NA = size(A,2); NZ = size(Z,2);
[tmp,T__on_ij,tmp] = intersect(transpose(1:MA),T__on{1+ncov},'stable');
[tmp,T_off_ij,tmp] = intersect(transpose(1:MA),T_off{1+ncov},'stable');
tmp__on__on = (A(T__on{1+ncov},:)*transpose(Z(S__on{1+ncov},:)));
ZRz__on__on = diag((tmp__on__on)*transpose(tmp__on__on)) - NA*ls__on;
ZCz__on__on = diag(transpose(A(T__on{1+ncov},:))*(tmp__on__on)*Z(S__on{1+ncov},:)) - lt__on*ls__on;
tmp__on_off = (A(T__on{1+ncov},:)*transpose(Z(S_off{1+ncov},:)));
ZRz__on_off = diag((tmp__on_off)*transpose(tmp__on_off)) - NA*ls_off;
ZCz__on_off = diag(transpose(A(T__on{1+ncov},:))*(tmp__on_off)*Z(S_off{1+ncov},:)) - lt__on*ls_off;
tmp_off__on = (A(T_off{1+ncov},:)*transpose(Z(S__on{1+ncov},:)));
ZRz_off__on = diag((tmp_off__on)*transpose(tmp_off__on)) - NA*ls__on;
ZCz_off__on = diag(transpose(A(T_off{1+ncov},:))*(tmp_off__on)*Z(S__on{1+ncov},:)) - lt_off*ls__on;
tmp_off_off = (A(T_off{1+ncov},:)*transpose(Z(S_off{1+ncov},:)));
ZRz_off_off = diag((tmp_off_off)*transpose(tmp_off_off)) - NA*ls_off;
ZCz_off_off = diag(transpose(A(T_off{1+ncov},:))*(tmp_off_off)*Z(S_off{1+ncov},:)) - lt_off*ls_off;
pz = (ls__on)/(ls__on+ls_off); qz = (ls_off)/(ls__on+ls_off);
ZRz__on = ...
        + min(ZRz__on__on*pz/pz,ZRz__on_off*pz/qz) ...
        + min(ZRz__on__on*qz/pz,ZRz__on_off*qz/qz) ...
  ;
ZRz_off = ...
        + min(ZRz_off__on*pz/pz,ZRz_off_off*pz/qz) ... 
        + min(ZRz_off__on*qz/pz,ZRz_off_off*qz/qz) ...
  ;
ZRz{1+ncov}(T__on_ij) = ZRz__on;
ZRz{1+ncov}(T_off_ij) = ZRz_off;
papz = lt__on/(lt__on+lt_off)*ls__on/(ls__on+ls_off); 
paqz = lt__on/(lt__on+lt_off)*ls_off/(ls__on+ls_off);
qapz = lt_off/(lt__on+lt_off)*ls__on/(ls__on+ls_off); 
qaqz = lt_off/(lt__on+lt_off)*ls_off/(ls__on+ls_off);
ZCz{1+ncov} = ...
  + min( ZCz__on__on*papz/papz , min( ZCz__on_off*papz/paqz , min( ZCz_off__on*papz/qapz , ZCz_off_off*papz/qaqz ))) ...
  + min( ZCz__on__on*paqz/papz , min( ZCz__on_off*paqz/paqz , min( ZCz_off__on*paqz/qapz , ZCz_off_off*paqz/qaqz ))) ...
  + min( ZCz__on__on*qapz/papz , min( ZCz__on_off*qapz/paqz , min( ZCz_off__on*qapz/qapz , ZCz_off_off*qapz/qaqz ))) ...
  + min( ZCz__on__on*qaqz/papz , min( ZCz__on_off*qaqz/paqz , min( ZCz_off__on*qaqz/qapz , ZCz_off_off*qaqz/qaqz ))) ...
  ;
ZRz{1+ncov} = ZRz{1+ncov} * (MA/MZ)*(NA/NZ)*(NA/NZ);
ZCz{1+ncov} = ZCz{1+ncov} * (MA/MZ)*(MA/MZ)*(NA/NZ);
end;%if ( cov_bother_flag_AAAA(1+ncov) &  cov_bother_flag_AZZA(1+ncov));

if (~cov_bother_flag_AAAA(1+ncov) &  cov_bother_flag_AZZA(1+ncov));
MA = lt__on + lt_off; MZ = ls__on + ls_off; NA = size(A,2); NZ = size(Z,2);
[tmp,T__on_ij,tmp] = intersect(transpose(1:MA),T__on{1+ncov},'stable');
[tmp,T_off_ij,tmp] = intersect(transpose(1:MA),T_off{1+ncov},'stable');
tmp__on__on = (A(T__on{1+ncov},:)*transpose(Z(S__on{1+ncov},:)));
ZRz__on__on = diag((tmp__on__on)*transpose(tmp__on__on)) - NA*ls__on;
ZCz__on__on = diag(transpose(A(T__on{1+ncov},:))*(tmp__on__on)*Z(S__on{1+ncov},:)) - lt__on*ls__on;
tmp__on_off = (A(T__on{1+ncov},:)*transpose(Z(S_off{1+ncov},:)));
ZRz__on_off = diag((tmp__on_off)*transpose(tmp__on_off)) - NA*ls_off;
ZCz__on_off = diag(transpose(A(T__on{1+ncov},:))*(tmp__on_off)*Z(S_off{1+ncov},:)) - lt__on*ls_off;
tmp_off__on = (A(T_off{1+ncov},:)*transpose(Z(S__on{1+ncov},:)));
ZRz_off__on = diag((tmp_off__on)*transpose(tmp_off__on)) - NA*ls__on;
ZCz_off__on = diag(transpose(A(T_off{1+ncov},:))*(tmp_off__on)*Z(S__on{1+ncov},:)) - lt_off*ls__on;
tmp_off_off = (A(T_off{1+ncov},:)*transpose(Z(S_off{1+ncov},:)));
ZRz_off_off = diag((tmp_off_off)*transpose(tmp_off_off)) - NA*ls_off;
ZCz_off_off = diag(transpose(A(T_off{1+ncov},:))*(tmp_off_off)*Z(S_off{1+ncov},:)) - lt_off*ls_off;
pz = (ls__on)/(ls__on+ls_off); qz = (ls_off)/(ls__on+ls_off);
ZRz__on = ...
        + min(ZRz__on__on*pz/pz,ZRz__on_off*pz/qz) ...
        + min(ZRz__on__on*qz/pz,ZRz__on_off*qz/qz) ...
  ;
ZRz_off = ...
        + min(ZRz_off__on*pz/pz,ZRz_off_off*pz/qz) ... 
        + min(ZRz_off__on*qz/pz,ZRz_off_off*qz/qz) ...
  ;
ZRz{1+ncov}(T__on_ij) = ZRz__on;
ZRz{1+ncov}(T_off_ij) = ZRz_off;
papz = lt__on/(lt__on+lt_off)*ls__on/(ls__on+ls_off); 
paqz = lt__on/(lt__on+lt_off)*ls_off/(ls__on+ls_off);
qapz = lt_off/(lt__on+lt_off)*ls__on/(ls__on+ls_off); 
qaqz = lt_off/(lt__on+lt_off)*ls_off/(ls__on+ls_off);
if (lt__on>0);
ZCz{1+ncov} = ...
  + min( ZCz__on__on*papz/papz , ZCz__on_off*papz/paqz ) ...
  + min( ZCz__on__on*paqz/papz , ZCz__on_off*paqz/paqz ) ...
  + min( ZCz__on__on*qapz/papz , ZCz__on_off*qapz/paqz ) ...
  + min( ZCz__on__on*qaqz/papz , ZCz__on_off*qaqz/paqz ) ...
  ;
end;%if (lt__on>0);
if (lt_off>0);
ZCz{1+ncov} = ...
  + min( ZCz_off__on*papz/qapz , ZCz_off_off*papz/qaqz ) ...
  + min( ZCz_off__on*paqz/qapz , ZCz_off_off*paqz/qaqz ) ...
  + min( ZCz_off__on*qapz/qapz , ZCz_off_off*qapz/qaqz ) ...
  + min( ZCz_off__on*qaqz/qapz , ZCz_off_off*qaqz/qaqz ) ...
  ;
end;%if (lt_off>0);
ZRz{1+ncov} = ZRz{1+ncov} * (MA/MZ)*(NA/NZ)*(NA/NZ);
ZCz{1+ncov} = ZCz{1+ncov} * (MA/MZ)*(MA/MZ)*(NA/NZ);
end;%if (~cov_bother_flag_AAAA(1+ncov) &  cov_bother_flag_AZZA(1+ncov));

if ( cov_bother_flag_AAAA(1+ncov) & ~cov_bother_flag_AZZA(1+ncov));
MA = lt__on + lt_off; MZ = ls__on + ls_off; NA = size(A,2); NZ = size(Z,2);
[tmp,T__on_ij,tmp] = intersect(transpose(1:MA),T__on{1+ncov},'stable');
[tmp,T_off_ij,tmp] = intersect(transpose(1:MA),T_off{1+ncov},'stable');
tmp__on__on = (A(T__on{1+ncov},:)*transpose(Z(S__on{1+ncov},:)));
ZRz__on__on = diag((tmp__on__on)*transpose(tmp__on__on)) - NA*ls__on;
ZCz__on__on = diag(transpose(A(T__on{1+ncov},:))*(tmp__on__on)*Z(S__on{1+ncov},:)) - lt__on*ls__on;
tmp__on_off = (A(T__on{1+ncov},:)*transpose(Z(S_off{1+ncov},:)));
ZRz__on_off = diag((tmp__on_off)*transpose(tmp__on_off)) - NA*ls_off;
ZCz__on_off = diag(transpose(A(T__on{1+ncov},:))*(tmp__on_off)*Z(S_off{1+ncov},:)) - lt__on*ls_off;
tmp_off__on = (A(T_off{1+ncov},:)*transpose(Z(S__on{1+ncov},:)));
ZRz_off__on = diag((tmp_off__on)*transpose(tmp_off__on)) - NA*ls__on;
ZCz_off__on = diag(transpose(A(T_off{1+ncov},:))*(tmp_off__on)*Z(S__on{1+ncov},:)) - lt_off*ls__on;
tmp_off_off = (A(T_off{1+ncov},:)*transpose(Z(S_off{1+ncov},:)));
ZRz_off_off = diag((tmp_off_off)*transpose(tmp_off_off)) - NA*ls_off;
ZCz_off_off = diag(transpose(A(T_off{1+ncov},:))*(tmp_off_off)*Z(S_off{1+ncov},:)) - lt_off*ls_off;
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
ZRz{1+ncov}(T__on_ij) = ZRz__on;
ZRz{1+ncov}(T_off_ij) = ZRz_off;
papz = lt__on/(lt__on+lt_off)*ls__on/(ls__on+ls_off); 
paqz = lt__on/(lt__on+lt_off)*ls_off/(ls__on+ls_off);
qapz = lt_off/(lt__on+lt_off)*ls__on/(ls__on+ls_off); 
qaqz = lt_off/(lt__on+lt_off)*ls_off/(ls__on+ls_off);
if (ls__on>0);
ZCz{1+ncov} = ...
  + min( ZCz__on__on*papz/papz , ZCz_off__on*papz/qapz ) ...
  + min( ZCz__on__on*paqz/papz , ZCz_off__on*paqz/qapz ) ...
  + min( ZCz__on__on*qapz/papz , ZCz_off__on*qapz/qapz ) ...
  + min( ZCz__on__on*qaqz/papz , ZCz_off__on*qaqz/qapz ) ...
  ;
end;%if (ls_on>0);
if (ls_off>0);
ZCz{1+ncov} = ...
  + min( ZCz__on_off*papz/paqz , ZCz_off_off*papz/qaqz ) ...
  + min( ZCz__on_off*paqz/paqz , ZCz_off_off*paqz/qaqz ) ...
  + min( ZCz__on_off*qapz/paqz , ZCz_off_off*qapz/qaqz ) ...
  + min( ZCz__on_off*qaqz/paqz , ZCz_off_off*qaqz/qaqz ) ...
  ;
end;%if (ls_off>0);
ZRz{1+ncov} = ZRz{1+ncov} * (MA/MZ)*(NA/NZ)*(NA/NZ);
ZCz{1+ncov} = ZCz{1+ncov} * (MA/MZ)*(MA/MZ)*(NA/NZ);
end;%if (~cov_bother_flag_AAAA(1+ncov) & cov_bother_flag_AZZA(1+ncov));

end;%for ncov=0:ncovs-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finished looping through covariates 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cov_set_AAAA = 0;
for ncov=0:ncovs-1;
if (cov_bother_flag_AAAA(1+ncov) & ~cov_set_AAAA);
ZRa_min = ZRa{1+ncov} ; ZCa_min = ZCa{1+ncov} ; cov_set_AAAA=1;
elseif (cov_bother_flag_AAAA(1+ncov) & cov_set_AAAA);
ZRa_min = min(ZRa_min,ZRa{1+ncov}) ; ZCa_min = min(ZCa_min,ZCa{1+ncov}) ; 
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
ZRz_min = ZRz{1+ncov}; ZCz_min = ZCz{1+ncov} ; cov_set_AZZA=1;
elseif ((cov_bother_flag_AAAA(1+ncov) | cov_bother_flag_AZZA(1+ncov)) & cov_set_AZZA);
ZRz_min = min(ZRz_min,ZRz{1+ncov}); ZCz_min = min(ZCz_min,ZCz{1+ncov}) ; 
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

ZR = ZRa_min(:) - ZRz_min(:) ; ZC = ZCa_min(:) - ZCz_min(:) ;
[tmp,tmp_rij] = sort(ZR,'ascend'); [tmp,tmp_cij] = sort(ZC,'ascend');
% if gamma=0, then we set up gamma_tmp_row to remove a single row ;
%if (gamma>0); gamma_tmp_row = max(gamma,(1-1e-6)/length(tmp_rij)); ammag_tmp_row = min(1-gamma,(length(tmp_rij)-1)/length(tmp_rij));
%elseif gamma<=0; gamma_tmp_row = (1-1e-6)/length(tmp_rij); ammag_tmp_row = (length(tmp_rij)-1)/length(tmp_rij);
%end;%if gamma==0;
% setting up ammag_tmp_col to remove as many cols as necessary so that log(ncols_pos)/log(nrows_pos) = log(ncols_pre)/log(nrows_pre) ;
% i.e., log(ammag_tmp_col*ncols_pre)/log(ammag_tmp_row*nrows_pre) = log(ncols_pre)/log(nrows_pre) ;
% i.e., log(ammag_tmp_col*ncols_pre) = (log(ammag_tmp_row) + log(nrows_pre))*log(ncols_pre)/log(nrows_pre) ;
% i.e., log(ammag_tmp_col) = log(ammag_tmp_row)*log(ncols_pre)/log(nrows_pre) ;
% i.e., ammag_tmp_col = exp(log(ammag_tmp_row)*log(ncols_pre)/log(nrows_pre));
%ammag_tmp_col = exp(log(ammag_tmp_row)*log(length(tmp_cij))/log(length(tmp_rij)));
%gamma_tmp_col = 1-ammag_tmp_col-1e-6;
%rdrop = r_rem{iteration}(tmp_rij(1:ceil(gamma_tmp_row*end))); cdrop = c_rem{iteration}(tmp_cij(1:ceil(gamma_tmp_col*end)));
%rkeep = r_rem{iteration}(tmp_rij(ceil(gamma_tmp_row*end)+1:end)); ckeep = c_rem{iteration}(tmp_cij(ceil(gamma_tmp_col*end)+1:end));
%r_rem{iteration+1} = rkeep;
%c_rem{iteration+1} = ckeep;
%rij=[rij,rdrop(end:-1:1)]; cij=[cij,cdrop(end:-1:1)];
[MA,NA] = size(A); [MZ,NZ] = size(Z);
tmp_R = sum(ZR); R_d = log(MA) + log(NA) + log(MA) + log(NA);
tmp_C = sum(ZC); C_d = log(NA) + log(MA) + log(NA) + log(MA);
out_trace(iteration,:) = [iteration , MA , NA , sign(tmp_R)*exp(log(abs(tmp_R))-R_d) , sign(tmp_C)*exp(log(abs(tmp_C))-C_d) sum(cov_bother_flag_AAAA)+sum(cov_bother_flag_AZZA)];
%for nr=1:length(rdrop);
%out_xdrop(out_xdrop_ij,:) = [A_n_rind_vals_lookup(rdrop(end+1-nr))-1 , -1];
%out_xdrop_ij = out_xdrop_ij+1;
%end;%for nr=1:length(rdrop);
%for nc=1:length(cdrop);
%out_xdrop(out_xdrop_ij,:) = [-1 , A_n_cind_vals_lookup(cdrop(end+1-nc))-1];
%out_xdrop_ij = out_xdrop_ij+1;
%end;%for nc=1:length(cdrop);
%A = A(tmp_rij(ceil(gamma_tmp_row*end)+1:end),tmp_cij(ceil(gamma_tmp_col*end)+1:end));
%Z = Z(:,tmp_cij(ceil(gamma_tmp_col*end)+1:end));
%T = T(tmp_rij(ceil(gamma_tmp_row*end)+1:end),:);
%S = S(:,:);
%iteration = iteration+1;
%end;%while min(size(A))>2;
%rij = [rij,setdiff(row_ij,rij)]; cij = [cij,setdiff(col_ij,cij)];
%rij = rij(end:-1:1); cij =  cij(end:-1:1);
LR = out_trace(1,4);
LC = out_trace(1,5);
Lx = out_trace(1,6);
Lx_AAAA = sum(cov_bother_flag_AAAA);
Lx_AZZA = sum(cov_bother_flag_AZZA);
