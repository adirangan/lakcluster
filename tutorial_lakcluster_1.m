function [rij,cij] = tutorial_lakcluster_1(fname_in,loopS_flag,loopS_nrmax,loopS_ncmax);
% Simple matlab implementation of categorical-covariate-corrected low-rank-biclustering routine. ;
% Instead of using low-rank updates, we simply recalculate the row- and col-scores at each iteration. ;
% Safe to use with one sided covariates; that is, covariates which extend across cases or controls, but not both ;
% Also, note that only T covariates are used (V ignored) ;
% loopS_flag==1 calculates the loopscore for each corner of the (1:nrmax,1:ncmax) submatrix of the sorted array ; 
% Note that this simple implementation assumes sparsity ~ 0.5 ;

if nargin<1;
disp('testing tutorial_lakcluster_1.m with tutorial_lakcluster_1_test.m');
tutorial_lakcluster_1_test;
return;
end;%if nargin<1;

if nargin<2; loopS_flag=0; loopS_nrmax=0; loopS_ncmax=0; end;
GLOBAL_A_n_repl_num=0;

fid = fopen(fname_in,'r');
continue_flag=1;
while continue_flag;
var = fscanf(fid,'%[^=]');
equ = fscanf(fid,'%c',1);
spc = fscanf(fid,'%c',1);
if strcmp(var,'GLOBAL_verbose'); GLOBAL_verbose = fscanf(fid,'%d'); end;
if strcmp(var,'GLOBAL_thread_count'); GLOBAL_thread_count = fscanf(fid,'%d'); end;
if strcmp(var,'GLOBAL_CFILTER_ERRCHECK'); GLOBAL_CFILTER_ERRCHECK = fscanf(fid,'%d'); end;
if strcmp(var,'GLOBAL_CFILTER_SCOREOUT'); GLOBAL_CFILTER_SCOREOUT = fscanf(fid,'%d'); end;
if strcmp(var,'GLOBAL_CFILTER_AGGRESSIVE'); GLOBAL_CFILTER_AGGRESSIVE = fscanf(fid,'%f'); end;
if strcmp(var,'GLOBAL_force_kr'); GLOBAL_force_kr = fscanf(fid,'%d'); end;
if strcmp(var,'GLOBAL_force_wk'); GLOBAL_force_wk = fscanf(fid,'%d'); end;
if strcmp(var,'GLOBAL_force_xv'); GLOBAL_force_xv = fscanf(fid,'%d'); end;
if strcmp(var,'A_n_name'); A_n_name = fscanf(fid,'%[^;]'); end;
if strcmp(var,'A_t_name'); A_t_name = fscanf(fid,'%[^;]'); end;
if strcmp(var,'A_n_rows'); A_n_rows = fscanf(fid,'%d'); end;
if strcmp(var,'A_n_cols'); A_n_cols = fscanf(fid,'%d'); end;
if strcmp(var,'A_n_rind'); A_n_rind = fscanf(fid,'%[^;]'); end;
if strcmp(var,'A_n_cind'); A_n_cind = fscanf(fid,'%[^;]'); end;
if strcmp(var,'Z_n_name'); Z_n_name = fscanf(fid,'%[^;]'); end;
if strcmp(var,'Z_t_name'); Z_t_name = fscanf(fid,'%[^;]'); end;
if strcmp(var,'Z_n_rows'); Z_n_rows = fscanf(fid,'%d'); end;
if strcmp(var,'Z_n_rind'); Z_n_rind = fscanf(fid,'%[^;]'); end;
if strcmp(var,'T_n_name'); T_n_name = fscanf(fid,'%[^;]'); end;
if strcmp(var,'T_t_name'); T_t_name = fscanf(fid,'%[^;]'); end;
if strcmp(var,'T_n_cols'); T_n_cols = fscanf(fid,'%d'); end;
if strcmp(var,'T_n_cind'); T_n_cind = fscanf(fid,'%[^;]'); end;
if strcmp(var,'S_n_name'); S_n_name = fscanf(fid,'%[^;]'); end;
if strcmp(var,'S_t_name'); S_t_name = fscanf(fid,'%[^;]'); end;
if strcmp(var,'V_n_name'); V_n_name = fscanf(fid,'%[^;]'); end;
if strcmp(var,'V_t_name'); V_t_name = fscanf(fid,'%[^;]'); end;
if strcmp(var,'V_n_rows'); V_n_rows = fscanf(fid,'%d'); end;
if strcmp(var,'V_n_rind'); V_n_rind = fscanf(fid,'%[^;]'); end;
if strcmp(var,'GLOBAL_prm_flag'); prm_flag = fscanf(fid,'%d'); end;
if strcmp(var,'GLOBAL_out_name'); GLOBAL_out_name = fscanf(fid,'%[^;]'); end;
if strcmp(var,'A_n_repl_num'); GLOBAL_A_n_repl_num = fscanf(fid,'%d'); end;
if strcmp(var,'A_n_repl_name');
GLOBAL_A_n_repl_name = cell(GLOBAL_A_n_repl_num,1);
for nl=1:GLOBAL_A_n_repl_num;
if (nl<GLOBAL_A_n_repl_num); GLOBAL_A_n_repl_name{nl} = fscanf(fid,'%[^,]'); smc = fscanf(fid,'%c',1); smc = fscanf(fid,'%c',1); elseif (nl==GLOBAL_A_n_repl_num); GLOBAL_A_n_repl_name{nl} = fscanf(fid,'%[^;]'); end;
end;%for nl=1:GLOBAL_A_n_repl_num;
end;%if strcmp(var,'A_n_repl_name');
if strcmp(var,'END'); continue_flag=0; end;
smc = fscanf(fid,'%c',1);
nln = fscanf(fid,'%c',1);
end;%while continue_flag;
fclose(fid);

A_n_rind_vals = tutorial_binary_uncompress(A_n_rind,1:A_n_rows,1)>0;
A_n_cind_vals = tutorial_binary_uncompress(A_n_cind,1:A_n_cols,1)>0;
Z_n_rind_vals = tutorial_binary_uncompress(Z_n_rind,1:Z_n_rows,1)>0;
T_n_cind_vals = tutorial_binary_uncompress(T_n_cind,1:T_n_cols,1)>0;
V_n_rind_vals = tutorial_binary_uncompress(V_n_rind,1:V_n_rows,1)>0;

A_n_rind_vals_lookup = find(A_n_rind_vals);
A_n_cind_vals_lookup = find(A_n_cind_vals);
Z_n_rind_vals_lookup = find(Z_n_rind_vals);
T_n_cind_vals_lookup = find(T_n_cind_vals);

A_orig = tutorial_binary_uncompress(A_n_name,find(A_n_rind_vals),find(A_n_cind_vals));
Z_orig = tutorial_binary_uncompress(Z_n_name,find(Z_n_rind_vals),find(A_n_cind_vals));
T_orig = tutorial_binary_uncompress(T_n_name,find(A_n_rind_vals),find(T_n_cind_vals))>0;
S_orig = tutorial_binary_uncompress(S_n_name,find(Z_n_rind_vals),find(T_n_cind_vals))>0;

if (GLOBAL_A_n_repl_num>0);
disp(sprintf(' %% found %d previous biclusters: ',GLOBAL_A_n_repl_num));
for nl=1:GLOBAL_A_n_repl_num; disp(sprintf(' %% %% %s ',GLOBAL_A_n_repl_name{nl})); end;%for nl=1:GLOBAL_A_n_repl_num;
for nl=1:GLOBAL_A_n_repl_num; 
bc_xdrop = textread(sprintf('%s',GLOBAL_A_n_repl_name{nl})); 
bc_rdrop = bc_xdrop(:,1); bc_cdrop = bc_xdrop(:,2);
[tmp_ij,bc_rdrop_ij,A_n_rind_ij] = intersect(1+bc_rdrop,find(A_n_rind_vals));
[tmp_ij,bc_cdrop_ij,A_n_cind_ij] = intersect(1+bc_cdrop,find(A_n_cind_vals));
A_orig(A_n_rind_ij,A_n_cind_ij) = 2*(randn(length(A_n_rind_ij),length(A_n_cind_ij))>0)-1;
[tmp_ij,bc_rdrop_ij,Z_n_rind_ij] = intersect(1+bc_rdrop,find(Z_n_rind_vals));
Z_orig(Z_n_rind_ij,A_n_cind_ij) = 2*(randn(length(Z_n_rind_ij),length(A_n_cind_ij))>0)-1;
disp(sprintf(' %% %% nl %d: replacing %d-x-%d in A and %d-x-%d in Z out of %d-x-%d',nl,length(A_n_rind_ij),length(A_n_cind_ij),length(Z_n_rind_ij),length(A_n_cind_ij),length(find(bc_rdrop>-1)),length(find(bc_cdrop>-1))));
end;%for nl=1:GLOBAL_A_n_repl_num;
end;%if (GLOBAL_A_n_repl_num>0);

if (loopS_flag>0);
Dinc_n_rij = find(A_n_rind_vals); 
Xinc_n_rij = find(Z_n_rind_vals);
g_ij = find(A_n_cind_vals);
read_out_xdrop = textread(sprintf('%s_xdrop.txt',GLOBAL_out_name));
rdrop = read_out_xdrop(:,1); cdrop = read_out_xdrop(:,2);
rij = rdrop(find(rdrop>=0)); rij = [rij(:) ; setdiff(Dinc_n_rij(:)-1,rij(:))]; rij = rij(end:-1:1);
cij = cdrop(find(cdrop>=0)); cij = [cij(:) ; setdiff(g_ij(:)-1,cij(:))]; cij = cij(end:-1:1);
[tmp,tmp,rij_lookdn] = intersect(rij,Dinc_n_rij(:)-1,'stable'); rij_lookdn=rij_lookdn-1; A_n_rind_vals_lookup = 1+rij; T_n_rind_vals_lookup = 1+rij;
[tmp,tmp,cij_lookdn] = intersect(cij,g_ij(:)-1,'stable'); cij_lookdn=cij_lookdn-1; A_n_cind_vals_lookup = 1+cij; Z_n_cind_vals_lookup = 1+cij;
A_orig = A_orig(1+rij_lookdn,1+cij_lookdn);
T_orig = T_orig(1+rij_lookdn,:);
Z_orig = Z_orig(:,1+cij_lookdn);
S_orig = S_orig(:,:);
A_orig = A_orig(1:loopS_nrmax,1:loopS_ncmax);
T_orig = T_orig(1:loopS_nrmax,:);
Z_orig = Z_orig(:,1:loopS_ncmax);
S_orig = S_orig(:,:);
end;%if (loopS_flag>0);

ncovs = size(T_orig,2)-1; % the first column of T,S should be all ones ;
gamma = GLOBAL_CFILTER_AGGRESSIVE;
if loopS_flag==0; disp(sprintf('beginning tutorial_lakcluster_1: gamma %f',gamma)); end;
if loopS_flag==1; disp(sprintf('beginning tutorial_lakcluster_1: loopS_nrmax %d loopS_ncmax %d',loopS_nrmax,loopS_ncmax)); end;
[MA,NA] = size(A_orig); [MZ,NZ] = size(Z_orig);
out_xdrop = zeros(MA+NA,2); out_xdrop_ij=1;
out_trace = zeros(MA+NA,6); 
row_ij = 1:MA; col_ij = 1:NA;
A = A_orig; Z = Z_orig; T = T_orig>0; S = S_orig>0;
rij = [];cij = [];
iteration=1;
r_rem{iteration} = row_ij;c_rem{iteration} = col_ij;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if loopS_flag==0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin iterations ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while min(size(A))>2;
disp(sprintf('iteration %d, A size %d-x-%d = exp(%0.2f)-x-exp(%0.2f) ratio %0.2f',iteration,size(A),log(size(A)),size(A,2)/size(A,1)));
ncovs = size(T,2)-1; % the first column of T,S should be all ones ;
cov_bother_flag_AAAA=zeros(1+ncovs,1);cov_bother_flag_AAAA(1)=1;
cov_bother_flag_AZZA=zeros(1+ncovs,1);cov_bother_flag_AZZA(1)=1;
clear ZRa ZCa ZRz ZCz ;
ZRa = cell(ncovs);ZCa = cell(ncovs);
ZRz = cell(ncovs);ZCz = cell(ncovs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through covariates ;
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
% Finished looping through covariates ;
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
if (gamma>0); gamma_tmp_row = max(gamma,(1-1e-6)/length(tmp_rij)); ammag_tmp_row = min(1-gamma,(length(tmp_rij)-1)/length(tmp_rij));
elseif gamma<=0; gamma_tmp_row = (1-1e-6)/length(tmp_rij); ammag_tmp_row = (length(tmp_rij)-1)/length(tmp_rij);
end;%if gamma==0;
% setting up ammag_tmp_col to remove as many cols as necessary so that log(ncols_pos)/log(nrows_pos) = log(ncols_pre)/log(nrows_pre) ;
% i.e., log(ammag_tmp_col*ncols_pre)/log(ammag_tmp_row*nrows_pre) = log(ncols_pre)/log(nrows_pre) ;
% i.e., log(ammag_tmp_col*ncols_pre) = (log(ammag_tmp_row) + log(nrows_pre))*log(ncols_pre)/log(nrows_pre) ;
% i.e., log(ammag_tmp_col) = log(ammag_tmp_row)*log(ncols_pre)/log(nrows_pre) ;
% i.e., ammag_tmp_col = exp(log(ammag_tmp_row)*log(ncols_pre)/log(nrows_pre));
ammag_tmp_col = exp(log(ammag_tmp_row)*log(length(tmp_cij))/log(length(tmp_rij)));
gamma_tmp_col = 1-ammag_tmp_col-1e-6;
rdrop = r_rem{iteration}(tmp_rij(1:ceil(gamma_tmp_row*end))); cdrop = c_rem{iteration}(tmp_cij(1:ceil(gamma_tmp_col*end)));
rkeep = r_rem{iteration}(tmp_rij(ceil(gamma_tmp_row*end)+1:end)); ckeep = c_rem{iteration}(tmp_cij(ceil(gamma_tmp_col*end)+1:end));
r_rem{iteration+1} = rkeep;
c_rem{iteration+1} = ckeep;
rij=[rij,rdrop(end:-1:1)]; cij=[cij,cdrop(end:-1:1)];
[MA,NA] = size(A); [MZ,NZ] = size(Z);
tmp_R = sum(ZR); R_d = log(MA) + log(NA) + log(MA) + log(NA);
tmp_C = sum(ZC); C_d = log(NA) + log(MA) + log(NA) + log(MA);
out_trace(iteration,:) = [iteration , MA , NA , sign(tmp_R)*exp(log(abs(tmp_R))-R_d) , sign(tmp_C)*exp(log(abs(tmp_C))-C_d) sum(cov_bother_flag_AAAA)+sum(cov_bother_flag_AZZA)];
for nr=1:length(rdrop);
out_xdrop(out_xdrop_ij,:) = [A_n_rind_vals_lookup(rdrop(end+1-nr))-1 , -1];
out_xdrop_ij = out_xdrop_ij+1;
end;%for nr=1:length(rdrop);
for nc=1:length(cdrop);
out_xdrop(out_xdrop_ij,:) = [-1 , A_n_cind_vals_lookup(cdrop(end+1-nc))-1];
out_xdrop_ij = out_xdrop_ij+1;
end;%for nc=1:length(cdrop);
A = A(tmp_rij(ceil(gamma_tmp_row*end)+1:end),tmp_cij(ceil(gamma_tmp_col*end)+1:end));
Z = Z(:,tmp_cij(ceil(gamma_tmp_col*end)+1:end));
T = T(tmp_rij(ceil(gamma_tmp_row*end)+1:end),:);
S = S(:,:);
iteration = iteration+1;
end;%while min(size(A))>2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End iterations ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rij = [rij,setdiff(row_ij,rij)]; cij = [cij,setdiff(col_ij,cij)];
rij = rij(end:-1:1); cij =  cij(end:-1:1);

fid = fopen(sprintf('%s_xdrop.txt',GLOBAL_out_name),'w');
fprintf(fid,'%d %d\n',transpose(out_xdrop(1:out_xdrop_ij-1,:)));
fclose(fid);
fid = fopen(sprintf('%s_trace.txt',GLOBAL_out_name),'w');
fprintf(fid,'%d %d %d %0.16f %0.16f %d\n',transpose(out_trace(1:iteration-1,:)));
fclose(fid);

plot_flag=0;
if plot_flag;
xpan = 128;
xpander = zeros(ncovs,ncovs*xpan); for nx=1:xpan; xpander(:,nx:xpan:end) = eye(ncovs); end;
figure;cla;
imagesc([ ...
A_orig(rij,cij) , zeros(length(rij),32) , repmat(T_orig(rij,2:end)*xpander,1,1) ; ...
zeros(32,length(cij)) , zeros(32,32) , zeros(32,ncovs*xpan) ; ...
Z_orig(:,cij) , zeros(size(Z_orig,1),32) , repmat(S_orig(:,2:end)*xpander,1,1) ; 
],[-1,1]*1.5);
title('reordered matrix');
end;%if plot_flag;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end;%if loopS_flag==0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if loopS_flag>0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nrmax = loopS_nrmax; ncmax = loopS_ncmax; nr_ra = 2:nrmax; nc_ra = floor(linspace(4,ncmax,min(256,ncmax))); nrbins = length(nr_ra); ncbins = length(nc_ra);
LR2_ra = zeros(nrbins,ncbins); LC2_ra = zeros(nrbins,ncbins); Lx2_ra = zeros(nrbins,ncbins);
%LR2_ra = zeros(size(A_orig)); LC2_ra = zeros(size(A_orig)); Lx2_ra = zeros(size(A_orig));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outer loop (row updates) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_est_start = tic;
for nr0=length(nr_ra):-1:1;
if (mod(nr0,16)==0); t_est_time_elapsed = floor(toc(t_est_start)); t_est_time_remaining = floor(t_est_time_elapsed*nr0/(size(A_orig,1)-nr0)); disp(sprintf(' %% nr0 %d: rows %.3d/%.3d, time_elapsed %ds, estimated_time_remaining %ds',nr0,nr_ra(nr0),size(A_orig,1),t_est_time_elapsed,t_est_time_remaining)); end;
nrows_rem = nr_ra(nr0);
A = A_orig(1:nr_ra(nr0),:); Z = Z_orig; T = T_orig(1:nr_ra(nr0),:)>0; S = S_orig>0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inner loop (col updates) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nc0=length(nc_ra):-1:1;
ncols_rem = nc_ra(nc0);
[LR,LC,Lx] = tutorial_lakcluster_0_getscore(A(:,1:nc_ra(nc0)),Z(:,1:nc_ra(nc0)),T(:,:),S(:,:));
LR2_ra(nr0,nc0) = LR; LC2_ra(nr0,nc0) = LC; Lx2_ra(nr0,nc0) = Lx;
% -------------------------------- ;
% finish inner loop ;
% -------------------------------- ;
end;%for nc0=length(nc_ra):-1:1;
% -------------------------------- ;
% finish outer loop ;
% -------------------------------- ;
end;%for nr0=length(nr_ra):-1:1;
LR_ra = LR2_ra(:,end);
LC_ra = LC2_ra(end,:);
Lx_ra = Lx2_ra(:,end);
save(sprintf('%s_loopS.mat',GLOBAL_out_name),'ncbins','nrbins','ncmax','nrmax','nc_ra','nr_ra','LR_ra','LC_ra','Lx_ra','LR2_ra','LC2_ra','Lx2_ra');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end;%if loopS_flag>0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

