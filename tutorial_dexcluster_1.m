function [rij,cij] = tutorial_dexcluster_1(fname_in,loopS_flag,loopS_nrmax,loopS_ncmax);
% searches for rank-0 biclusters using auc. Designed for the no-covariate case (although will produce output when used with covariates). ;
% Assumes continuous (distinct) values; does not correct for nonunique values (although will produce output when given nonunique data points). ;
% When loopS_flag is set to 1 this calculates the loopscore for each corner of the (1:nrmax,1:ncmax) submatrix of the sorted array ; 

if nargin<1;
disp('testing tutorial_dexcluster_1.m with tutorial_dexcluster_1_test.m');
tutorial_dexcluster_1_test;
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
if loopS_flag==0; disp(sprintf('beginning tutorial_dexcluster_1: gamma %f',gamma)); end;
if loopS_flag==1; disp(sprintf('beginning tutorial_dexcluster_1: loopS_nrmax %d loopS_ncmax %d',loopS_nrmax,loopS_ncmax)); end;
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
std_A = transpose(unique(cov_cat_A)); std_Z = transpose(unique(cov_cat_Z));
if     (length(std_A)==1 & length(std_Z)==1); std_A = std_A; std_Z = std_Z; NSTD = 1;
elseif (length(std_A)> 1 & length(std_Z)==1); std_A = std_A; std_Z = repmat(std_Z,length(std_A),1); NSTD = length(std_A);
elseif (length(std_A)==1 & length(std_Z)> 1); std_A = repmat(std_A,length(std_Z),1); std_Z = std_Z; NSTD = length(std_Z);
elseif (length(std_A)> 1 & length(std_Z)> 1); std_A = union(std_A,std_Z); std_Z = std_A; NSTD = length(std_A);
end;% if;
std_A_ij_orig=cell(NSTD,1); std_Z_ij_orig=cell(NSTD,1);
for nstd=1:NSTD;
std_A_ij_orig{nstd} = find(cov_cat_A==std_A(nstd)); std_Z_ij_orig{nstd} = find(cov_cat_Z==std_Z(nstd));
end;%for nstd=1:NSTD;
Zcs = cell(NSTD,1);
Zir = cell(NSTD,1);
Auc_T = zeros(NA,NSTD);
Auc_P = cell(NSTD,1);
for nstd=1:NSTD;
Zcs{nstd} = zeros(length(std_A_ij_orig{nstd}) + length(std_Z_ij_orig{nstd}),NA);
Zir{nstd} = zeros(length(std_A_ij_orig{nstd}) + length(std_Z_ij_orig{nstd}),NA);
Auc_P{nstd} = 0.5 + zeros(length(std_A_ij_orig{nstd}),NA);
for nc=1:NA;
tmp_vec = [A(std_A_ij_orig{nstd},nc) ; Z(std_Z_ij_orig{nstd},nc)];
tmp_Aij = zeros(size(tmp_vec)); tmp_Aij(1:length(std_A_ij_orig{nstd}))=1; 
tmp_Zij = zeros(size(tmp_vec)); tmp_Zij(length(std_A_ij_orig{nstd}) + (1:length(std_Z_ij_orig{nstd})))=1;
[tmp_val,tmp_ij] = sort(tmp_vec,'ascend'); [tmp_val,tmp_ir] = sort(tmp_ij,'ascend');
tmp_Aij = tmp_Aij(tmp_ij); tmp_Zij = tmp_Zij(tmp_ij);
tmp_Zcs = cumsum(tmp_Zij);
tmp_auc = (mean(find(tmp_Aij))-mean(find(tmp_Zij)))/length(tmp_vec) + 0.5;
% if either std_A_ij_orig{nstd} or std_Z_ij_orig{nstd} is empty, we need to fix tmp_auc ;
if ~isfinite(tmp_auc); tmp_auc = 0.5; end;
% at this point length(find(Z(std_Z_ij_orig{nstd},nc)<A(std_A_ij_orig{nstd}(ij2),nc))) should equal tmp_Zcs(tmp_ir(ij2)) for all ij2;
% test with: 
% ij2=1;disp(sprintf(' %% %d vs %d',length(find(Z(std_Z_ij_orig{nstd},nc)<A(std_A_ij_orig{nstd}(ij2),nc))),tmp_Zcs(tmp_ir(ij2))));
Zcs{nstd}(:,nc) = tmp_Zcs;
Zir{nstd}(:,nc) = tmp_ir;
Auc_T(nc,nstd) = tmp_auc;
for nl=1:length(std_A_ij_orig{nstd}); Auc_P{nstd}(nl,nc) = tmp_Zcs(tmp_ir(nl))/length(std_Z_ij_orig{nstd}); if ~isfinite(Auc_P{nstd}(nl,nc)); Auc_P{nstd}(nl,nc) = 0.5; end; end;
if test_flag;
tmp_norm = 0; 
for nl=1:length(std_A_ij_orig{nstd}); 
tmp_norm = tmp_norm + (length(find(Z(std_Z_ij_orig{nstd},nc)<A(std_A_ij_orig{nstd}(nl),nc))) - tmp_Zcs(tmp_ir(nl))).^2; 
end;%for nl=1:length(std_A_ij_orig{nstd}); 
disp(sprintf(' %% cs-vs-ir error: %f',tmp_norm));
end;%if test_flag;
end;%for nc=1:NA;
end;%for nstd=1:NSTD;
if test_flag;
tmp_norm = 0; 
for nstd=1:NSTD; for nc=1:NA; for nl=1:length(std_A_ij_orig{nstd}); 
tmp_norm = tmp_norm + (length(find(Z(std_Z_ij_orig{nstd},nc)<A(std_A_ij_orig{nstd}(nl),nc))) - Zcs{nstd}(Zir{nstd}(nl,nc),nc)).^2; 
end;end;end;%for nstd=1:NSTD; for nc=1:NA; for nl=1:length(std_A_ij_orig{nstd}); 
disp(sprintf(' %% cs-vs-ir error: %f',tmp_norm));
end;%if test_flag;
Auc_P_cmb = zeros(MA,1);
fill_flag = ones(MA,1);
for nstd=1:NSTD;
%Auc_P_cmb(std_A_ij_orig{nstd}) = sum(2*abs(Auc_P{nstd}(:,:)-0.5),2);
tmp_row = std_A_ij_orig{nstd};
tmp_sum = sum(2*abs(Auc_P{nstd}(:,:)-0.5),2);
[tmp_rij,tmp_rij_a,tmp_rij_b] = intersect(tmp_row,find(fill_flag),'stable');
Auc_P_cmb(tmp_rij) = tmp_sum(tmp_rij_a);
fill_flag(tmp_rij)=0;
end;%for nstd=1:NSTD;
Auc_T_orig = Auc_T; Auc_P_cmb_orig = Auc_P_cmb; Auc_P_orig = Auc_P;
std_A_ij = std_A_ij_orig;
std_Z_ij = std_Z_ij_orig;
ncovs_rem = 0; for nstd=1:NSTD; ncovs_rem = ncovs_rem + (length(std_A_ij{nstd})>0); end;
nstd=1;nrows_rem = length(std_A_ij{nstd}); for nstd=2:NSTD; nrows_rem = min(nrows_rem,length(std_A_ij{nstd})); end;
nstd=1;nrows_sum = length(std_A_ij{nstd}); for nstd=2:NSTD; nrows_sum = nrows_sum + length(std_A_ij{nstd}); end;
ncols_rem = NA;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if loopS_flag==0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin iterations ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while min(nrows_sum,ncols_rem)>2;
disp(sprintf('iteration %d (ncovs_rem %d), A size %d-x-%d = exp(%0.2f)-x-exp(%0.2f) ratio %0.2f',iteration,ncovs_rem,nrows_sum,ncols_rem,log(nrows_sum),log(ncols_rem),ncols_rem/nrows_sum));
ZC = 2*min(abs(Auc_T(c_rem{iteration},:)-0.5),[],2);
ZR = Auc_P_cmb(r_rem{iteration});
if test_flag;
ZC2 = zeros(size(A,2),1); 
for nc2=1:size(A,2);
[tmp_avg,tmp_min,tmp_auc] = getcauc_2(A(:,nc2),Z(:,nc2),[],[],T*transpose(2.^(0:ncovs)),S*transpose(2.^(0:ncovs)));
ZC2(nc2) = 2*abs(tmp_min-0.5);
end;%for nc2=1:size(A,2);
ZR2 = zeros(size(A,1),size(A,2));
for nr2=1:size(A,1);
for nc2=1:size(A,2);
tmp_cov_A = T*transpose(2.^(0:ncovs)); tmp_cov_Z = S*transpose(2.^(0:ncovs));
tmp_ij = find(tmp_cov_Z==tmp_cov_A(nr2));
ZR2(nr2,nc2) = length(find(A(nr2,nc2)>Z(tmp_ij,nc2)))/length(tmp_ij);
end;%for nc2=1:size(A,2);
end;%for nr2=1:size(A,1);
ZR2 = sum(2*abs(ZR2-0.5),2);
subplot(1,2,1);cla;plot(1:length(ZC),ZC,'b.-',1:length(ZC2),ZC2,'bo-');
title(sprintf('col score: iter %d',iteration));
subplot(1,2,2);cla;plot(1:length(ZR),ZR,'r.-',1:length(ZR2),ZR2,'ro-');
title(sprintf('row score: iter %d',iteration));
drawnow; pause();
end;%if test_flag;
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
% for each dropped row we modify Auc_T ;
for nr=1:length(rdrop);
nstd_ra = find(std_A == cov_cat_A(rdrop(nr)));
for nstd_ij = 1:length(nstd_ra);
tmp_list = std_A_ij_orig{nstd_ra(nstd_ij)};
tmp_ij = find(tmp_list==rdrop(nr));
for nc=1:NA;
tmp_list = std_A_ij{nstd_ra(nstd_ij)};
Auc_T(nc,nstd_ra(nstd_ij)) = (Auc_T(nc,nstd_ra(nstd_ij))*(length(tmp_list)) - Auc_P{nstd_ra(nstd_ij)}(tmp_ij,nc))/(length(tmp_list)-1);
end;%for nc=1:NA;
tmp_list = std_A_ij{nstd_ra(nstd_ij)};
std_A_ij{nstd_ra(nstd_ij)} = setdiff(tmp_list,rdrop(nr));
end;%for nstd_ij = 1:length(nstd_ra);
end;%for nr=1:length(rdrop);
% for each dropped col we modify Auc_P_cmb ;
for nc=1:length(cdrop);
for nr=1:MA;
nstd_ra = find(std_A == cov_cat_A(nr),1,'first');
for nstd_ij = 1:length(nstd_ra);
tmp_list = std_A_ij_orig{nstd_ra(nstd_ij)};
tmp_ij = find(tmp_list==nr);
Auc_P_cmb(nr) = Auc_P_cmb(nr) - 2*abs(Auc_P{nstd_ra(nstd_ij)}(tmp_ij,cdrop(nc))-0.5);
end;%for nstd_ij = 1:length(nstd_ra);
end;%for nr=1:MA;
end;%for nc=1:length(cdrop);
r_rem{iteration+1} = rkeep;
c_rem{iteration+1} = ckeep;
rij=[rij,rdrop(end:-1:1)]; cij=[cij,cdrop(end:-1:1)];
tmp_R = mean(ZR)/ncols_rem;
tmp_C = mean(ZC);
out_trace(iteration,:) = [iteration , nrows_sum , ncols_rem , tmp_R , tmp_C , ncovs_rem];
for nr=1:length(rdrop);
out_xdrop(out_xdrop_ij,:) = [A_n_rind_vals_lookup(rdrop(end+1-nr))-1 , -1];
out_xdrop_ij = out_xdrop_ij+1;
end;%for nr=1:length(rdrop);
for nc=1:length(cdrop);
out_xdrop(out_xdrop_ij,:) = [-1 , A_n_cind_vals_lookup(cdrop(end+1-nc))-1];
out_xdrop_ij = out_xdrop_ij+1;
end;%for nc=1:length(cdrop);
ncovs_rem = 0; for nstd=1:NSTD; tmp_list = std_A_ij{nstd}; ncovs_rem = ncovs_rem + (length(tmp_list)>0); end;
nstd=1;tmp_list = std_A_ij{nstd};nrows_rem = length(tmp_list); for nstd=2:NSTD; tmp_list = std_A_ij{nstd}; nrows_rem = min(nrows_rem,length(tmp_list)); end;
nstd=1;tmp_list = std_A_ij{nstd};nrows_sum = length(tmp_list); for nstd=2:NSTD; tmp_list = std_A_ij{nstd}; nrows_sum = nrows_sum + length(tmp_list); end;
ncols_rem = ncols_rem - length(cdrop);
if test_flag;
A = A(tmp_rij(ceil(gamma_tmp_row*end)+1:end),tmp_cij(ceil(gamma_tmp_col*end)+1:end));
Z = Z(:,tmp_cij(ceil(gamma_tmp_col*end)+1:end));
T = T(tmp_rij(ceil(gamma_tmp_row*end)+1:end),:);
S = S(:,:);
end;%if test_flag;
iteration = iteration+1;
end;%while min(nrows_sum,ncols_rem)>2;

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

nrmax = loopS_nrmax; ncmax = loopS_ncmax;
nr_ra = 1:nrmax; nc_ra = 1:ncmax;
Auc_T = Auc_T_orig; Auc_P_cmb = Auc_P_cmb_orig; Auc_P = Auc_P_orig;
LR2_ra = zeros(size(A_orig)); LC2_ra = zeros(size(A_orig)); Lx2_ra = zeros(size(A_orig));
t_est_start = tic;
for nr0=size(A_orig,1):-1:1;
if (mod(nr0,16)==0); t_est_time_elapsed = floor(toc(t_est_start)); t_est_time_remaining = floor(t_est_time_elapsed*nr0/(size(A_orig,1)-nr0)); disp(sprintf(' %% nr0 %d: rows %.3d/%.3d, time_elapsed %ds, estimated_time_remaining %ds',nr0,nr_ra(nr0),size(A_orig,1),t_est_time_elapsed,t_est_time_remaining)); end;
Auc_T_tmp = Auc_T; Auc_P_cmb_tmp = Auc_P_cmb; Auc_P_tmp = Auc_P;
ncovs_rem_tmp = 0; for nstd=1:NSTD; tmp_list = std_A_ij{nstd}; ncovs_rem_tmp = ncovs_rem_tmp + (length(tmp_list)>0); end;
for nc0=size(A_orig,2):-1:1;
ZC_tmp = 2*min(abs(Auc_T_tmp(1:nc0,:)-0.5),[],2);
LC2_ra(nr0,nc0) = mean(ZC_tmp);
ZR_tmp = Auc_P_cmb_tmp(1:nr0);
LR2_ra(nr0,nc0) = mean(ZR_tmp)/nc0;
Lx2_ra(nr0,nc0) = ncovs_rem_tmp;
% for each dropped col we modify Auc_P_cmb_tmp ;
cdrop = nc0;
for nr=1:nr0;
nstd_ra = find(std_A == cov_cat_A(nr),1,'first');
for nstd_ij = 1:length(nstd_ra);
tmp_list = std_A_ij_orig{nstd_ra(nstd_ij)};
tmp_ij = find(tmp_list==nr);
Auc_P_cmb_tmp(nr) = Auc_P_cmb_tmp(nr) - 2*abs(Auc_P_tmp{nstd_ra(nstd_ij)}(tmp_ij,cdrop)-0.5);
end;%for nstd_ij = 1:length(nstd_ra);
end;%for nr=1:nr0;
end;%for nc0=size(A_orig,2):-1:1;
% for each dropped row we modify Auc_T ;
rdrop = nr0;
nstd_ra = find(std_A == cov_cat_A(rdrop));
for nstd_ij = 1:length(nstd_ra);
tmp_list = std_A_ij_orig{nstd_ra(nstd_ij)};
tmp_ij = find(tmp_list==rdrop);
for nc=1:NA;
tmp_list = std_A_ij{nstd_ra(nstd_ij)};
Auc_T(nc,nstd_ra(nstd_ij)) = (Auc_T(nc,nstd_ra(nstd_ij))*(length(tmp_list)) - Auc_P{nstd_ra(nstd_ij)}(tmp_ij,nc))/(length(tmp_list)-1);
end;%for nc=1:NA;
tmp_list = std_A_ij{nstd_ra(nstd_ij)};
std_A_ij{nstd_ra(nstd_ij)} = setdiff(tmp_list,rdrop);
end;%for nstd_ij = 1:length(nstd_ra);
end;%for nr0=size(A_orig,1):-1:1;

nrbins = size(A_orig,1); ncbins = size(A_orig,2);
LR_ra = LR2_ra(:,end);
LC_ra = LC2_ra(end,:);
Lx_ra = Lx2_ra(:,end);
save(sprintf('%s_loopS.mat',GLOBAL_out_name),'ncbins','nrbins','ncmax','nrmax','nc_ra','nr_ra','LR_ra','LC_ra','Lx_ra','LR2_ra','LC2_ra','Lx2_ra');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end;%if loopS_flag>0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
