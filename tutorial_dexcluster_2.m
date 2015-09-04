function [rij,cij] = tutorial_dexcluster_2(fname_in,loopS_flag,loopS_nrmax,loopS_ncmax);
% searches for rank-0 biclusters using a modified version of covariate-corrected auc (cauc). ;
% Covariate-correction is implemented. ;
% We also allow one sided covariates; that is, covariates which extend across cases or controls, but not both ;
% Assumes continuous (distinct) values; does not correct for nonunique values (although will produce output when given nonunique data points). ;
% Also, note that only T covariates are used (V ignored);
% When loopS_flag is set to 1 this calculates the loopscore for each corner of the (1:nrmax,1:ncmax) submatrix of the sorted array ; 

if nargin<1;
disp('testing tutorial_dexcluster_2.m with tutorial_dexcluster_2__test.m');
tutorial_dexcluster_2_test;
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

load(A_n_name,'MATRIX_ORIG');
A_orig = MATRIX_ORIG(find(A_n_rind_vals),find(A_n_cind_vals));
Z_orig = MATRIX_ORIG(find(Z_n_rind_vals),find(A_n_cind_vals));
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
for nl=1:length(A_n_cind_ij);
tmp_Z = Z_orig(:,A_n_cind_ij(nl)); tmp_A = interp1(1:length(tmp_Z),sort(tmp_Z,'ascend'),1+(length(tmp_Z)-1)*rand(length(A_n_rind_ij),1));
A_orig(A_n_rind_ij,A_n_cind_ij(nl)) = tmp_A;
end;%for nl=1:length(A_n_cind_ij);
[tmp_ij,bc_rdrop_ij,Z_n_rind_ij] = intersect(1+bc_rdrop,find(Z_n_rind_vals));
disp(sprintf(' %% %% nl %d: replacing %d-x-%d in A and leaving %d-x-%d in Z (out of %d-x-%d)',nl,length(A_n_rind_ij),length(A_n_cind_ij),length(Z_n_rind_ij),length(A_n_cind_ij),length(find(bc_rdrop>-1)),length(find(bc_cdrop>-1))));
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

test_flag=0;
% Test starting here with: ;
% MA=132; NA=531; MZ = floor(MA/2)+1; A_orig=randn(MA,NA); A_orig(1:40,1:60) = A_orig(1:40,1:60)+5.5; Z_orig=randn(MZ,NA); T_orig=ones(size(A_orig,1),2); S_orig=ones(size(Z_orig,1),2); T_orig(:,2)=randn(size(A_orig,1),1)>0; S_orig(:,2)=randn(size(Z_orig,1),1)>0; GLOBAL_CFILTER_AGGRESSIVE=0; A_n_rind_vals_lookup = 1:MA; A_n_cind_vals_lookup = 1:NA; test_flag=0;
ncovs = size(T_orig,2)-1; % the first column of T,S should be all ones ;
gamma = GLOBAL_CFILTER_AGGRESSIVE;
if loopS_flag==0; disp(sprintf('beginning tutorial_dexcluster_2: gamma %f',gamma)); end;
if loopS_flag==1; disp(sprintf('beginning tutorial_dexcluster_2: loopS_nrmax %d loopS_ncmax %d',loopS_nrmax,loopS_ncmax)); end;
[MA,NA] = size(A_orig); [MZ,NZ] = size(Z_orig);
out_xdrop = zeros(MA+NA,2); out_xdrop_ij=1;
out_trace = zeros(MA+NA,6); 
row_ij = 1:MA; col_ij = 1:NA;
A = A_orig; Z = Z_orig; T = T_orig>0; S = S_orig>0;
rij = [];cij = [];
iteration=1;
r_rem{iteration} = row_ij;c_rem{iteration} = col_ij;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preliminary setup ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[MA,NA] = size(A); [MZ,NZ] = size(Z); 
ncovs = size(T,2)-1; % the first column of T,S should be all ones ;
cov_cat_A = T*transpose(2.^(0:ncovs)); cov_cat_Z = S*transpose(2.^(0:ncovs));
cov_A_up_orig_=cell(ncovs,1); 
cov_A_dn_orig_=cell(ncovs,1); 
cov_Z_up_orig_=cell(ncovs,1);
cov_Z_dn_orig_=cell(ncovs,1);
for ncov=0:ncovs-1;
cov_A_up_orig_{1+ncov} = find(T(:,1+1+ncov)==1);
cov_A_dn_orig_{1+ncov} = find(T(:,1+1+ncov)==0);
cov_Z_up_orig_{1+ncov} = find(S(:,1+1+ncov)==1);
cov_Z_dn_orig_{1+ncov} = find(S(:,1+1+ncov)==0);
end;%for ncov=0:ncovs-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Zcs_all = zeros(MA+MZ,NA); Zir_all = zeros(MA+MZ,NA); X_P_all = 0.5 + zeros(MA,NA); X_T_all = zeros(NA,1);
for nc=1:NA;
tmp_vec = [A(:,nc) ; Z(:,nc)]; tmp_Aij = zeros(size(tmp_vec)); tmp_Aij(1:MA)=1; tmp_Zij = zeros(size(tmp_vec)); tmp_Zij(MA + (1:MZ))=1;
[tmp_val,tmp_ij] = sort(tmp_vec,'ascend'); [tmp_val,tmp_ir] = sort(tmp_ij,'ascend');
tmp_Aij = tmp_Aij(tmp_ij); tmp_Zij = tmp_Zij(tmp_ij); tmp_Zcs = cumsum(tmp_Zij);
tmp_auc = (mean(find(tmp_Aij))-mean(find(tmp_Zij)))/length(tmp_vec) + 0.5;
% if either block of tmp_vec is empty, we need to fix tmp_auc ;
if ~isfinite(tmp_auc); tmp_auc = 0.5; end;
% test with: tmp_x=0; for nra=1:MA;for nrz=1:MZ; tmp_x=tmp_x+(A(nra,nc)>Z(nrz,nc)); end;end;
% tmp_X = (mean(find(tmp_Aij))-mean(find(tmp_Zij)) + 0.5*(MA+MZ))*(MA*MZ)/(MA+MZ);
% tmp_X = ((sum(find(tmp_Aij))*MZ-sum(find(tmp_Zij)*MA))/(MA+MZ) + 0.5*MA*MZ);
tmp_X = sum(tmp_Zcs(find(tmp_Aij)));
% at this point length(find(Z(:,nc)<A(ij2,nc))) should equal tmp_Zcs(tmp_ir(ij2)) for all ij2;
% test with: for ij2=1:MA;disp(sprintf(' %% %d vs %d',length(find(Z(:,nc)<A(ij2,nc))),tmp_Zcs(tmp_ir(ij2))));end;
Zcs_all(:,nc) = tmp_Zcs;
Zir_all(:,nc) = tmp_ir;
X_T_all(nc) = tmp_X;
% Note that, with this construction, X_T_all = sum(X_P_all,1);
for nl=1:MA; X_P_all(nl,nc) = tmp_Zcs(tmp_ir(nl)); if ~isfinite(X_P_all(nl,nc)); X_P_all(nl,nc) = 0; end; end;
%if test_flag; tmp_norm = 0; for nl=1:MA; tmp_norm = tmp_norm + (length(find(Z(:,nc)<A(nl,nc))) - tmp_Zcs(tmp_ir(nl))).^2; end; disp(sprintf(' %% cs-vs-ir error: %f',tmp_norm)); end;%if test_flag;
end;%for nc=1:NA;
if test_flag; tmp_norm = 0; for nc=1:NA; for nl=1:MA; tmp_norm = tmp_norm + (length(find(Z(:,nc)<A(nl,nc))) - X_P_all(nl,nc)).^2; end; end; disp(sprintf(' %% _all cs-vs-ir error: %f',tmp_norm)); end;%if test_flag;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Zcs_xx_up = cell(ncovs,1); Zir_xx_up = cell(ncovs,1); X_T_xx_up = zeros(NA,ncovs); X_P_xx_up = cell(ncovs,1);
Zcs_xx_dn = cell(ncovs,1); Zir_xx_dn = cell(ncovs,1); X_T_xx_dn = zeros(NA,ncovs); X_P_xx_dn = cell(ncovs,1);
Zcs_up_up = cell(ncovs,1); Zir_up_up = cell(ncovs,1); X_T_up_up = zeros(NA,ncovs); X_P_up_up = cell(ncovs,1);
Zcs_up_dn = cell(ncovs,1); Zir_up_dn = cell(ncovs,1); X_T_up_dn = zeros(NA,ncovs); X_P_up_dn = cell(ncovs,1);
Zcs_dn_up = cell(ncovs,1); Zir_dn_up = cell(ncovs,1); X_T_dn_up = zeros(NA,ncovs); X_P_dn_up = cell(ncovs,1);
Zcs_dn_dn = cell(ncovs,1); Zir_dn_dn = cell(ncovs,1); X_T_dn_dn = zeros(NA,ncovs); X_P_dn_dn = cell(ncovs,1);
for ncov=1:ncovs;
Zcs_xx_up{ncov} = zeros(MA + length(cov_Z_up_orig_{ncov}),NA); Zir_xx_up{ncov} = zeros(MA + length(cov_Z_up_orig_{ncov}),NA); X_P_xx_up{ncov} = zeros(MA,NA); 
Zcs_xx_dn{ncov} = zeros(MA + length(cov_Z_dn_orig_{ncov}),NA); Zir_xx_dn{ncov} = zeros(MA + length(cov_Z_dn_orig_{ncov}),NA); X_P_xx_dn{ncov} = zeros(MA,NA); 
Zcs_up_up{ncov} = zeros(length(cov_A_up_orig_{ncov}) + length(cov_Z_up_orig_{ncov}),NA); Zir_up_up{ncov} = zeros(length(cov_A_up_orig_{ncov}) + length(cov_Z_up_orig_{ncov}),NA); X_P_up_up{ncov} = zeros(length(cov_A_up_orig_{ncov}),NA); 
Zcs_up_dn{ncov} = zeros(length(cov_A_up_orig_{ncov}) + length(cov_Z_dn_orig_{ncov}),NA); Zir_up_dn{ncov} = zeros(length(cov_A_up_orig_{ncov}) + length(cov_Z_dn_orig_{ncov}),NA); X_P_up_dn{ncov} = zeros(length(cov_A_up_orig_{ncov}),NA); 
Zcs_dn_up{ncov} = zeros(length(cov_A_dn_orig_{ncov}) + length(cov_Z_up_orig_{ncov}),NA); Zir_dn_up{ncov} = zeros(length(cov_A_dn_orig_{ncov}) + length(cov_Z_up_orig_{ncov}),NA); X_P_dn_up{ncov} = zeros(length(cov_A_dn_orig_{ncov}),NA); 
Zcs_dn_dn{ncov} = zeros(length(cov_A_dn_orig_{ncov}) + length(cov_Z_dn_orig_{ncov}),NA); Zir_dn_dn{ncov} = zeros(length(cov_A_dn_orig_{ncov}) + length(cov_Z_dn_orig_{ncov}),NA); X_P_dn_dn{ncov} = zeros(length(cov_A_dn_orig_{ncov}),NA); 
for nc=1:NA;
tmp_vec = [A(1:MA,nc) ; Z(cov_Z_up_orig_{ncov},nc)]; tmp_Aij = zeros(size(tmp_vec)); tmp_Aij(1:MA)=1; tmp_Zij = zeros(size(tmp_vec)); tmp_Zij(MA + (1:length(cov_Z_up_orig_{ncov})))=1;
[tmp_val,tmp_ij] = sort(tmp_vec,'ascend'); [tmp_val,tmp_ir] = sort(tmp_ij,'ascend');
tmp_Aij = tmp_Aij(tmp_ij); tmp_Zij = tmp_Zij(tmp_ij); tmp_Zcs = cumsum(tmp_Zij);
tmp_auc = (mean(find(tmp_Aij))-mean(find(tmp_Zij)))/length(tmp_vec) + 0.5; if ~isfinite(tmp_auc); tmp_auc = 0.5; end;
tmp_X = sum(tmp_Zcs(find(tmp_Aij)));
Zcs_xx_up{ncov}(:,nc) = tmp_Zcs; Zir_xx_up{ncov}(:,nc) = tmp_ir; X_T_xx_up(nc,ncov) = tmp_X;
for nl=1:MA; X_P_xx_up{ncov}(nl,nc) = tmp_Zcs(tmp_ir(nl)); if ~isfinite(X_P_xx_up{ncov}(nl,nc)); X_P_xx_up{ncov}(nl,nc) = 0; end; end;
%if test_flag; tmp_norm = 0; for nl=1:MA; tmp_norm = tmp_norm + (length(find(Z(cov_Z_up_orig_{ncov},nc)<A(nl,nc))) - tmp_Zcs(tmp_ir(nl))).^2; end; disp(sprintf(' %% cs-vs-ir error: %f',tmp_norm)); end;%if test_flag;
tmp_vec = [A(1:MA,nc) ; Z(cov_Z_dn_orig_{ncov},nc)]; tmp_Aij = zeros(size(tmp_vec)); tmp_Aij(1:MA)=1; tmp_Zij = zeros(size(tmp_vec)); tmp_Zij(MA + (1:length(cov_Z_dn_orig_{ncov})))=1;
[tmp_val,tmp_ij] = sort(tmp_vec,'ascend'); [tmp_val,tmp_ir] = sort(tmp_ij,'ascend');
tmp_Aij = tmp_Aij(tmp_ij); tmp_Zij = tmp_Zij(tmp_ij); tmp_Zcs = cumsum(tmp_Zij);
tmp_auc = (mean(find(tmp_Aij))-mean(find(tmp_Zij)))/length(tmp_vec) + 0.5; if ~isfinite(tmp_auc); tmp_auc = 0.5; end;
tmp_X = ((sum(find(tmp_Aij))*length(cov_Z_dn_orig_{ncov})-sum(find(tmp_Zij)*MA))/(MA+length(cov_Z_dn_orig_{ncov})) + 0.5*MA*length(cov_Z_dn_orig_{ncov}));
Zcs_xx_dn{ncov}(:,nc) = tmp_Zcs; Zir_xx_dn{ncov}(:,nc) = tmp_ir; X_T_xx_dn(nc,ncov) = tmp_X;
for nl=1:MA; X_P_xx_dn{ncov}(nl,nc) = tmp_Zcs(tmp_ir(nl)); if ~isfinite(X_P_xx_dn{ncov}(nl,nc)); X_P_xx_dn{ncov}(nl,nc) = 0; end; end;
%if test_flag; tmp_norm = 0; for nl=1:MA; tmp_norm = tmp_norm + (length(find(Z(cov_Z_dn_orig_{ncov},nc)<A(nl,nc))) - tmp_Zcs(tmp_ir(nl))).^2; end; disp(sprintf(' %% cs-vs-ir error: %f',tmp_norm)); end;%if test_flag;
tmp_vec = [A(cov_A_up_orig_{ncov},nc) ; Z(cov_Z_up_orig_{ncov},nc)]; 
tmp_Aij = zeros(size(tmp_vec)); tmp_Aij(1:length(cov_A_up_orig_{ncov}))=1; tmp_Zij = zeros(size(tmp_vec)); tmp_Zij(length(cov_A_up_orig_{ncov}) + (1:length(cov_Z_up_orig_{ncov})))=1;
[tmp_val,tmp_ij] = sort(tmp_vec,'ascend'); [tmp_val,tmp_ir] = sort(tmp_ij,'ascend');
tmp_Aij = tmp_Aij(tmp_ij); tmp_Zij = tmp_Zij(tmp_ij); tmp_Zcs = cumsum(tmp_Zij);
tmp_auc = (mean(find(tmp_Aij))-mean(find(tmp_Zij)))/length(tmp_vec) + 0.5; if ~isfinite(tmp_auc); tmp_auc = 0.5; end;
tmp_X = sum(tmp_Zcs(find(tmp_Aij)));
Zcs_up_up{ncov}(:,nc) = tmp_Zcs; Zir_up_up{ncov}(:,nc) = tmp_ir; X_T_up_up(nc,ncov) = tmp_X;
for nl=1:length(cov_A_up_orig_{ncov}); X_P_up_up{ncov}(nl,nc) = tmp_Zcs(tmp_ir(nl)); if ~isfinite(X_P_up_up{ncov}(nl,nc)); X_P_up_up{ncov}(nl,nc) = 0; end; end;
%if test_flag; tmp_norm = 0; for nl=1:length(cov_A_up_orig_{ncov}); tmp_norm = tmp_norm + (length(find(Z(cov_Z_up_orig_{ncov},nc)<A(cov_A_up_orig_{ncov}(nl),nc))) - tmp_Zcs(tmp_ir(nl))).^2; end; disp(sprintf(' %% cs-vs-ir error: %f',tmp_norm)); end;%if test_flag;
tmp_vec = [A(cov_A_up_orig_{ncov},nc) ; Z(cov_Z_dn_orig_{ncov},nc)]; 
tmp_Aij = zeros(size(tmp_vec)); tmp_Aij(1:length(cov_A_up_orig_{ncov}))=1; tmp_Zij = zeros(size(tmp_vec)); tmp_Zij(length(cov_A_up_orig_{ncov}) + (1:length(cov_Z_dn_orig_{ncov})))=1;
[tmp_val,tmp_ij] = sort(tmp_vec,'ascend'); [tmp_val,tmp_ir] = sort(tmp_ij,'ascend');
tmp_Aij = tmp_Aij(tmp_ij); tmp_Zij = tmp_Zij(tmp_ij); tmp_Zcs = cumsum(tmp_Zij);
tmp_auc = (mean(find(tmp_Aij))-mean(find(tmp_Zij)))/length(tmp_vec) + 0.5; if ~isfinite(tmp_auc); tmp_auc = 0.5; end;
tmp_X = sum(tmp_Zcs(find(tmp_Aij)));
Zcs_up_dn{ncov}(:,nc) = tmp_Zcs; Zir_up_dn{ncov}(:,nc) = tmp_ir; X_T_up_dn(nc,ncov) = tmp_X;
for nl=1:length(cov_A_up_orig_{ncov}); X_P_up_dn{ncov}(nl,nc) = tmp_Zcs(tmp_ir(nl)); if ~isfinite(X_P_up_dn{ncov}(nl,nc)); X_P_up_dn{ncov}(nl,nc) = 0; end; end;
%if test_flag; tmp_norm = 0; for nl=1:length(cov_A_up_orig_{ncov}); tmp_norm = tmp_norm + (length(find(Z(cov_Z_dn_orig_{ncov},nc)<A(cov_A_up_orig_{ncov}(nl),nc))) - tmp_Zcs(tmp_ir(nl))).^2; end; disp(sprintf(' %% cs-vs-ir error: %f',tmp_norm)); end;%if test_flag;
tmp_vec = [A(cov_A_dn_orig_{ncov},nc) ; Z(cov_Z_up_orig_{ncov},nc)]; 
tmp_Aij = zeros(size(tmp_vec)); tmp_Aij(1:length(cov_A_dn_orig_{ncov}))=1; tmp_Zij = zeros(size(tmp_vec)); tmp_Zij(length(cov_A_dn_orig_{ncov}) + (1:length(cov_Z_up_orig_{ncov})))=1;
[tmp_val,tmp_ij] = sort(tmp_vec,'ascend'); [tmp_val,tmp_ir] = sort(tmp_ij,'ascend');
tmp_Aij = tmp_Aij(tmp_ij); tmp_Zij = tmp_Zij(tmp_ij); tmp_Zcs = cumsum(tmp_Zij);
tmp_auc = (mean(find(tmp_Aij))-mean(find(tmp_Zij)))/length(tmp_vec) + 0.5; if ~isfinite(tmp_auc); tmp_auc = 0.5; end;
tmp_X = sum(tmp_Zcs(find(tmp_Aij)));
Zcs_dn_up{ncov}(:,nc) = tmp_Zcs; Zir_dn_up{ncov}(:,nc) = tmp_ir; X_T_dn_up(nc,ncov) = tmp_X;
for nl=1:length(cov_A_dn_orig_{ncov}); X_P_dn_up{ncov}(nl,nc) = tmp_Zcs(tmp_ir(nl)); if ~isfinite(X_P_dn_up{ncov}(nl,nc)); X_P_dn_up{ncov}(nl,nc) = 0; end; end;
%if test_flag; tmp_norm = 0; for nl=1:length(cov_A_dn_orig_{ncov}); tmp_norm = tmp_norm + (length(find(Z(cov_Z_up_orig_{ncov},nc)<A(cov_A_dn_orig_{ncov}(nl),nc))) - tmp_Zcs(tmp_ir(nl))).^2; end; disp(sprintf(' %% cs-vs-ir error: %f',tmp_norm)); end;%if test_flag;
tmp_vec = [A(cov_A_dn_orig_{ncov},nc) ; Z(cov_Z_dn_orig_{ncov},nc)]; 
tmp_Aij = zeros(size(tmp_vec)); tmp_Aij(1:length(cov_A_dn_orig_{ncov}))=1; tmp_Zij = zeros(size(tmp_vec)); tmp_Zij(length(cov_A_dn_orig_{ncov}) + (1:length(cov_Z_dn_orig_{ncov})))=1;
[tmp_val,tmp_ij] = sort(tmp_vec,'ascend'); [tmp_val,tmp_ir] = sort(tmp_ij,'ascend');
tmp_Aij = tmp_Aij(tmp_ij); tmp_Zij = tmp_Zij(tmp_ij); tmp_Zcs = cumsum(tmp_Zij);
tmp_auc = (mean(find(tmp_Aij))-mean(find(tmp_Zij)))/length(tmp_vec) + 0.5; if ~isfinite(tmp_auc); tmp_auc = 0.5; end;
tmp_X = sum(tmp_Zcs(find(tmp_Aij)));
Zcs_dn_dn{ncov}(:,nc) = tmp_Zcs; Zir_dn_dn{ncov}(:,nc) = tmp_ir; X_T_dn_dn(nc,ncov) = tmp_X;
for nl=1:length(cov_A_dn_orig_{ncov}); X_P_dn_dn{ncov}(nl,nc) = tmp_Zcs(tmp_ir(nl)); if ~isfinite(X_P_dn_dn{ncov}(nl,nc)); X_P_dn_dn{ncov}(nl,nc) = 0; end; end;
%if test_flag; tmp_norm = 0; for nl=1:length(cov_A_dn_orig_{ncov}); tmp_norm = tmp_norm + (length(find(Z(cov_Z_dn_orig_{ncov},nc)<A(cov_A_dn_orig_{ncov}(nl),nc))) - tmp_Zcs(tmp_ir(nl))).^2; end; disp(sprintf(' %% cs-vs-ir error: %f',tmp_norm)); end;%if test_flag;
end;%for nc=1:NA;
end;%for ncov=1:ncovs;
if test_flag; tmp_norm = 0; 
for ncov=1:ncovs; for nc=1:NA; for nl=1:MA; tmp_norm = tmp_norm + (length(find(Z(cov_Z_up_orig_{ncov},nc)<A(nl,nc))) - Zcs_xx_up{ncov}(Zir_xx_up{ncov}(nl,nc),nc)).^2; end;end;end;%for ncov=1:ncovs; for nc=1:NA; for nl=1:MA; 
for ncov=1:ncovs; for nc=1:NA; for nl=1:MA; tmp_norm = tmp_norm + (length(find(Z(cov_Z_dn_orig_{ncov},nc)<A(nl,nc))) - Zcs_xx_dn{ncov}(Zir_xx_dn{ncov}(nl,nc),nc)).^2; end;end;end;%for ncov=1:ncovs; for nc=1:NA; for nl=1:MA; 
for ncov=1:ncovs; for nc=1:NA; for nl=1:length(cov_A_up_orig_{ncov}); tmp_norm = tmp_norm + (length(find(Z(cov_Z_up_orig_{ncov},nc)<A(cov_A_up_orig_{ncov}(nl),nc))) - Zcs_up_up{ncov}(Zir_up_up{ncov}(nl,nc),nc)).^2; end;end;end;%for ncov=1:ncovs; for nc=1:NA; for nl=1:length(cov_A_up_orig_{ncov}); 
for ncov=1:ncovs; for nc=1:NA; for nl=1:length(cov_A_up_orig_{ncov}); tmp_norm = tmp_norm + (length(find(Z(cov_Z_dn_orig_{ncov},nc)<A(cov_A_up_orig_{ncov}(nl),nc))) - Zcs_up_dn{ncov}(Zir_up_dn{ncov}(nl,nc),nc)).^2; end;end;end;%for ncov=1:ncovs; for nc=1:NA; for nl=1:length(cov_A_up_orig_{ncov}); 
for ncov=1:ncovs; for nc=1:NA; for nl=1:length(cov_A_dn_orig_{ncov}); tmp_norm = tmp_norm + (length(find(Z(cov_Z_up_orig_{ncov},nc)<A(cov_A_dn_orig_{ncov}(nl),nc))) - Zcs_dn_up{ncov}(Zir_dn_up{ncov}(nl,nc),nc)).^2; end;end;end;%for ncov=1:ncovs; for nc=1:NA; for nl=1:length(cov_A_dn_orig_{ncov}); 
for ncov=1:ncovs; for nc=1:NA; for nl=1:length(cov_A_dn_orig_{ncov}); tmp_norm = tmp_norm + (length(find(Z(cov_Z_dn_orig_{ncov},nc)<A(cov_A_dn_orig_{ncov}(nl),nc))) - Zcs_dn_dn{ncov}(Zir_dn_dn{ncov}(nl,nc),nc)).^2; end;end;end;%for ncov=1:ncovs; for nc=1:NA; for nl=1:length(cov_A_dn_orig_{ncov}); 
disp(sprintf(' %% ncov cs-vs-ir error: %f',tmp_norm));
end;%if test_flag;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear Zcs_xx_up; clear  Zir_xx_up;
clear Zcs_xx_dn; clear  Zir_xx_dn;
clear Zcs_up_up; clear  Zir_up_up;
clear Zcs_up_dn; clear  Zir_up_dn;
clear Zcs_dn_up; clear  Zir_dn_up;
clear Zcs_dn_dn; clear  Zir_dn_dn;
cov_A_up_=cov_A_up_orig_;
cov_A_dn_=cov_A_dn_orig_;
cov_Z_up_=cov_Z_up_orig_;
cov_Z_dn_=cov_Z_dn_orig_;
nrows_rem = MA;
ncols_rem = NA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if loopS_flag==0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin iterations ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while min(nrows_rem,ncols_rem)>2;
disp(sprintf('iteration %d, A size %d-x-%d = exp(%0.2f)-x-exp(%0.2f) ratio %0.2f',iteration,nrows_rem,ncols_rem,log(nrows_rem),log(ncols_rem),ncols_rem/nrows_rem));

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
% Here we assume that the X_P_* will be updated by modifying existing rows/cols to be filled with 0 ;
% This means that Z?z_g_* will be calculated correctly by summing, but Z?z_l will need care ;
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
ZR = sum(ZRz_min(r_rem{iteration},c_rem{iteration}),2); ZR = ZR(:); 
ZC = ZCz_min(:); ZC = ZC(c_rem{iteration});

if test_flag; [MA,NA] = size(A); tmp_norm = 0; for nc=1:NA; for nl=1:MA; tmp_norm = tmp_norm + (length(find(Z(:,nc)<A(nl,nc))) - X_P_all(r_rem{iteration}(nl),c_rem{iteration}(nc))).^2; end; end; disp(sprintf(' %% _all cs-vs-ir error: %f',tmp_norm)); end;%if test_flag; tmp_norm = 0; 
if test_flag; [MA,NA] = size(A); tmp_norm = 0; 
for ncov=1:ncovs; for nc=1:NA; for nl=1:MA; tmp_norm = tmp_norm + (length(find(Z(cov_Z_up_orig_{ncov},nc)<A(nl,nc))) - X_P_xx_up{ncov}(r_rem{iteration}(nl),c_rem{iteration}(nc))).^2; end;end;end;%for ncov=1:ncovs; for nc=1:NA; for nl=1:MA; 
for ncov=1:ncovs; for nc=1:NA; for nl=1:MA; tmp_norm = tmp_norm + (length(find(Z(cov_Z_dn_orig_{ncov},nc)<A(nl,nc))) - X_P_xx_dn{ncov}(r_rem{iteration}(nl),c_rem{iteration}(nc))).^2; end;end;end;%for ncov=1:ncovs; for nc=1:NA; for nl=1:MA; 
for ncov=1:ncovs; [tmp,cov_A_up_ij,tmp] = intersect(cov_A_up_orig_{ncov},cov_A_up_{ncov},'stable'); for nc=1:NA; for nl=1:length(cov_A_up_{ncov}); tmp_norm = tmp_norm + (length(find(Z(cov_Z_up_orig_{ncov},nc)<A_orig(cov_A_up_{ncov}(nl),c_rem{iteration}(nc)))) - X_P_up_up{ncov}(cov_A_up_ij(nl),c_rem{iteration}(nc))).^2; end;end;end;%for ncov=1:ncovs; for nc=1:NA; for nl=1:length(cov_A_up_{ncov}); 
for ncov=1:ncovs; [tmp,cov_A_up_ij,tmp] = intersect(cov_A_up_orig_{ncov},cov_A_up_{ncov},'stable'); for nc=1:NA; for nl=1:length(cov_A_up_{ncov}); tmp_norm = tmp_norm + (length(find(Z(cov_Z_dn_orig_{ncov},nc)<A_orig(cov_A_up_{ncov}(nl),c_rem{iteration}(nc)))) - X_P_up_dn{ncov}(cov_A_up_ij(nl),c_rem{iteration}(nc))).^2; end;end;end;%for ncov=1:ncovs; for nc=1:NA; for nl=1:length(cov_A_up_{ncov}); 
for ncov=1:ncovs; [tmp,cov_A_dn_ij,tmp] = intersect(cov_A_dn_orig_{ncov},cov_A_dn_{ncov},'stable'); for nc=1:NA; for nl=1:length(cov_A_dn_{ncov}); tmp_norm = tmp_norm + (length(find(Z(cov_Z_up_orig_{ncov},nc)<A_orig(cov_A_dn_{ncov}(nl),c_rem{iteration}(nc)))) - X_P_dn_up{ncov}(cov_A_dn_ij(nl),c_rem{iteration}(nc))).^2; end;end;end;%for ncov=1:ncovs; for nc=1:NA; for nl=1:length(cov_A_dn_{ncov}); 
for ncov=1:ncovs; [tmp,cov_A_dn_ij,tmp] = intersect(cov_A_dn_orig_{ncov},cov_A_dn_{ncov},'stable'); for nc=1:NA; for nl=1:length(cov_A_dn_{ncov}); tmp_norm = tmp_norm + (length(find(Z(cov_Z_dn_orig_{ncov},nc)<A_orig(cov_A_dn_{ncov}(nl),c_rem{iteration}(nc)))) - X_P_dn_dn{ncov}(cov_A_dn_ij(nl),c_rem{iteration}(nc))).^2; end;end;end;%for ncov=1:ncovs; for nc=1:NA; for nl=1:length(cov_A_dn_{ncov}); 
disp(sprintf(' %% ncov cs-vs-ir error: %f',tmp_norm));
end;%if test_flag;

% once ZR,ZC are determined ;
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
% for each dropped row we modify X_T_all, and for each ncov we modify both X_T_xx_up(nc,ncov) and X_T_xx_dn(nc,ncov), and some of the X_T_up_up, X_T_up_dn, X_T_dn_up, X_T_dn_dn ;
for nr=1:length(rdrop);
X_T_all(:) = X_T_all(:) - transpose(X_P_all(rdrop(nr),:)); X_P_all(rdrop(nr),:)=0;
for ncov=1:ncovs;
X_T_xx_up(:,ncov) = X_T_xx_up(:,ncov) - transpose(X_P_xx_up{ncov}(rdrop(nr),:)); X_P_xx_up{ncov}(rdrop(nr),:)=0;
X_T_xx_dn(:,ncov) = X_T_xx_dn(:,ncov) - transpose(X_P_xx_dn{ncov}(rdrop(nr),:)); X_P_xx_dn{ncov}(rdrop(nr),:)=0;
tmp_ij = find(cov_A_up_orig_{ncov}==rdrop(nr));
if ~isempty(tmp_ij); 
X_T_up_up(:,ncov) = X_T_up_up(:,ncov) - transpose(X_P_up_up{ncov}(tmp_ij,:)); X_P_up_up{ncov}(tmp_ij,:)=0;
X_T_up_dn(:,ncov) = X_T_up_dn(:,ncov) - transpose(X_P_up_dn{ncov}(tmp_ij,:)); X_P_up_dn{ncov}(tmp_ij,:)=0;
cov_A_up_{ncov} = setdiff(cov_A_up_{ncov},rdrop(nr));
end;% if found;
tmp_ij = find(cov_A_dn_orig_{ncov}==rdrop(nr));
if ~isempty(tmp_ij); 
X_T_dn_up(:,ncov) = X_T_dn_up(:,ncov) - transpose(X_P_dn_up{ncov}(tmp_ij,:)); X_P_dn_up{ncov}(tmp_ij,:)=0;
X_T_dn_dn(:,ncov) = X_T_dn_dn(:,ncov) - transpose(X_P_dn_dn{ncov}(tmp_ij,:)); X_P_dn_dn{ncov}(tmp_ij,:)=0;
cov_A_dn_{ncov} = setdiff(cov_A_dn_{ncov},rdrop(nr));
end;% if found;
end;%for ncov=1:ncovs;
end;%for nr=1:length(rdrop);
% for each dropped col we modify X_P_all, and for each ncov we modify X_P_xx_up, X_P_xx_dn, and some of the X_P_up_up, X_P_up_dn, X_P_dn_up, X_P_dn_dn ;
for nc=1:length(cdrop);
X_P_all(:,cdrop(nc)) = 0;
for ncov=1:ncovs;
X_P_xx_up{ncov}(:,cdrop(nc)) = 0;
X_P_xx_dn{ncov}(:,cdrop(nc)) = 0;
X_P_up_up{ncov}(:,cdrop(nc)) = 0;
X_P_up_dn{ncov}(:,cdrop(nc)) = 0;
X_P_dn_up{ncov}(:,cdrop(nc)) = 0;
X_P_dn_dn{ncov}(:,cdrop(nc)) = 0;
end;%for ncov=1:ncovs;
end;%for nc=1:length(cdrop);
r_rem{iteration+1} = rkeep;
c_rem{iteration+1} = ckeep;
rij=[rij,rdrop(end:-1:1)]; cij=[cij,cdrop(end:-1:1)];
MA = nrows_rem; NA = ncols_rem; MZ = size(Z_orig,1); NZ = ncols_rem;
tmp_R = sum(ZR); R_d = log(MA) + log(MZ) + log(NA);
tmp_C = sum(ZC); C_d = log(MA) + log(MZ) + log(NA);
tmp_R = sign(tmp_R)*exp(log(abs(tmp_R))-R_d);
tmp_C = sign(tmp_C)*exp(log(abs(tmp_C))-C_d);
out_trace(iteration,:) = [iteration , nrows_rem , ncols_rem , tmp_R , tmp_C , sum(cov_bother_flag_AAAA)+sum(cov_bother_flag_AZZA)];
for nr=1:length(rdrop);
out_xdrop(out_xdrop_ij,:) = [A_n_rind_vals_lookup(rdrop(end+1-nr))-1 , -1];
out_xdrop_ij = out_xdrop_ij+1;
end;%for nr=1:length(rdrop);
for nc=1:length(cdrop);
out_xdrop(out_xdrop_ij,:) = [-1 , A_n_cind_vals_lookup(cdrop(end+1-nc))-1];
out_xdrop_ij = out_xdrop_ij+1;
end;%for nc=1:length(cdrop);
ncols_rem = ncols_rem - length(cdrop);
nrows_rem = nrows_rem - length(rdrop);
if test_flag;
A = A(tmp_rij(ceil(gamma_tmp_row*end)+1:end),tmp_cij(ceil(gamma_tmp_col*end)+1:end));
Z = Z(:,tmp_cij(ceil(gamma_tmp_col*end)+1:end));
T = T(tmp_rij(ceil(gamma_tmp_row*end)+1:end),:);
S = S(:,:);
end;%if test_flag;
iteration = iteration+1;

end;%while min(nrows_rem,ncols_rem)>2;

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
tutorial_dexcluster_excerpt_0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outer loop (row updates) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_est_start = tic;
for nr0=length(nr_ra):-1:1;
if (mod(nr0,16)==0); t_est_time_elapsed = floor(toc(t_est_start)); t_est_time_remaining = floor(t_est_time_elapsed*nr0/(size(A_orig,1)-nr0)); disp(sprintf(' %% nr0 %d: rows %.3d/%.3d, time_elapsed %ds, estimated_time_remaining %ds',nr0,nr_ra(nr0),size(A_orig,1),t_est_time_elapsed,t_est_time_remaining)); end;
nrows_rem = nr_ra(nr0);
A = A_orig(1:nr_ra(nr0),:); Z = Z_orig; T = T_orig(1:nr_ra(nr0),:)>0; S = S_orig>0;
% -------------------------------- ;
% store backup ;
% -------------------------------- ;
cov_A_up_step_ = cov_A_up_;cov_A_dn_step_ = cov_A_dn_;cov_Z_up_step_ = cov_Z_up_;cov_Z_dn_step_ = cov_Z_dn_;
X_P_all_step = X_P_all;X_P_xx_up_step = X_P_xx_up;X_P_xx_dn_step = X_P_xx_dn;X_P_up_up_step = X_P_up_up;X_P_up_dn_step = X_P_up_dn;X_P_dn_up_step = X_P_dn_up;X_P_dn_dn_step = X_P_dn_dn;X_T_all_step = X_T_all;X_T_xx_up_step = X_T_xx_up;X_T_xx_dn_step = X_T_xx_dn;X_T_up_up_step = X_T_up_up;X_T_up_dn_step = X_T_up_dn;X_T_dn_up_step = X_T_dn_up;X_T_dn_dn_step = X_T_dn_dn;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inner loop (col updates) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nc0=length(nc_ra):-1:1;
ncols_rem = nc_ra(nc0);
% -------------------------------- ;
% calculate scores ;
% -------------------------------- ;
tutorial_dexcluster_excerpt_2;
tutorial_dexcluster_excerpt_t;
% -------------------------------- ;
% define drop & update auxiliary ;
% -------------------------------- ;
if nc0>1; MA = nrows_rem; NA = ncols_rem; tmp_rij_drop_full = []; tmp_rij_keep_full = [1:MA]; tmp_cij_drop_full = [nc_ra(nc0-1)+1:nc_ra(nc0)]; tmp_cij_keep_full = [1:nc_ra(nc0-1)]; tutorial_dexcluster_excerpt_3; end;% if nc0>1;
% -------------------------------- ;
% update LR2,LC2,Lx2 ;
% -------------------------------- ;
MA = nrows_rem; NA = ncols_rem; MZ = size(Z_orig,1); NZ = ncols_rem;
tmp_R = sum(ZR); R_d = log(MA) + log(MZ) + log(NA);
tmp_C = sum(ZC); C_d = log(MA) + log(MZ) + log(NA);
LR2_ra(nr0,nc0) = sign(tmp_R)*exp(log(abs(tmp_R))-R_d);
LC2_ra(nr0,nc0) = sign(tmp_C)*exp(log(abs(tmp_C))-C_d);
Lx2_ra(nr0,nc0) = sum(cov_bother_flag_AAAA)+sum(cov_bother_flag_AZZA);
% -------------------------------- ;
% update A,Z,T,S ;
% -------------------------------- ;
tutorial_dexcluster_excerpt_4;
% -------------------------------- ;
% finish inner loop ;
% -------------------------------- ;
end;%for nc0=length(nc_ra):-1:1;
% -------------------------------- ;
% reload backup ;
% -------------------------------- ;
A = A_orig(1:nr_ra(nr0),:); Z = Z_orig; T = T_orig(1:nr_ra(nr0),:)>0; S = S_orig>0;
cov_A_up_ = cov_A_up_step_;cov_A_dn_ = cov_A_dn_step_;cov_Z_up_ = cov_Z_up_step_;cov_Z_dn_ = cov_Z_dn_step_;
X_P_all = X_P_all_step;X_P_xx_up = X_P_xx_up_step;X_P_xx_dn = X_P_xx_dn_step;X_P_up_up = X_P_up_up_step;X_P_up_dn = X_P_up_dn_step;X_P_dn_up = X_P_dn_up_step;X_P_dn_dn = X_P_dn_dn_step;X_T_all = X_T_all_step;X_T_xx_up = X_T_xx_up_step;X_T_xx_dn = X_T_xx_dn_step;X_T_up_up = X_T_up_up_step;X_T_up_dn = X_T_up_dn_step;X_T_dn_up = X_T_dn_up_step;X_T_dn_dn = X_T_dn_dn_step;
% -------------------------------- ;
% update auxiliary data ;
% -------------------------------- ;
if nr0>1; MA = nrows_rem; NA = size(A_orig,2); tmp_rij_drop_full = [nr_ra(nr0-1)+1:nr_ra(nr0)]; tmp_rij_keep_full = [1:nr_ra(nr0-1)]; tmp_cij_drop_full = []; tmp_cij_keep_full = [1:NA]; tutorial_dexcluster_excerpt_3; end;% if nr0>1;
% -------------------------------- ;
% update A,Z,T,S ;
% -------------------------------- ;
tutorial_dexcluster_excerpt_4;
% -------------------------------- ;
% finish outer loop ;
% -------------------------------- ;
end;%for nr0=length(nr_ra):-1:1;
LR_ra = LR2_ra(:,end);
LC_ra = LC2_ra(end,:);
Lx_ra = Lx2_ra(:,end);
LR2_ra = 2*(LR2_ra-0.5); LC2_ra = 2*(LC2_ra-0.5);
LR_ra = 2*(LR_ra-0.5); LC_ra = 2*(LC_ra-0.5);

save(sprintf('%s_loopS.mat',GLOBAL_out_name),'ncbins','nrbins','ncmax','nrmax','nc_ra','nr_ra','LR_ra','LC_ra','Lx_ra','LR2_ra','LC2_ra','Lx2_ra');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end;%if loopS_flag>0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
