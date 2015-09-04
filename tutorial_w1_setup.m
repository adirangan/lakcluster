%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up a few flags and load data within tutorial_w1.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp_flag=0; clear_covariates_flag=0; load_B_flag=0; lowrank_vs_diffexp_flag=0;
path_gen = sprintf('%s%s.mat',path_base,gen_fname); path_gen2 = []; 
b_n0_ = strfind(path_gen,'_n0_'); if ~isempty(b_n0_); path_gen2 = path_gen; path_gen2(b_n0_+2) = '0'; path_gen2(b_n0_+3) = 'x'; clear_covariates_flag=1; lowrank_vs_diffexp_flag=+1; end;
b_n0x = strfind(path_gen,'_n0x'); if ~isempty(b_n0x); path_gen2 = path_gen; path_gen2(b_n0x+2) = '0'; lowrank_vs_diffexp_flag=+1; end;
b_n1_ = strfind(path_gen,'_n1_'); if ~isempty(b_n1_); path_gen2 = path_gen; path_gen2(b_n1_+2) = '1'; path_gen2(b_n1_+3) = 'x'; clear_covariates_flag=1; lowrank_vs_diffexp_flag=+1; end;
b_n1x = strfind(path_gen,'_n1x'); if ~isempty(b_n1x); path_gen2 = path_gen; path_gen2(b_n1x+2) = '1'; lowrank_vs_diffexp_flag=+1; end;
b_n2_ = strfind(path_gen,'_n2_'); if ~isempty(b_n2_); path_gen2 = path_gen; path_gen2(b_n2_+2) = '2'; path_gen2(b_n2_+3) = 'x'; clear_covariates_flag=1; load_B_flag=1; lowrank_vs_diffexp_flag=+1; end;
b_n2x = strfind(path_gen,'_n2x'); if ~isempty(b_n2x); path_gen2 = path_gen; path_gen2(b_n2x+2) = '2'; load_B_flag=1; lowrank_vs_diffexp_flag=+1; end;
b_nA_ = strfind(path_gen,'_nA_'); if ~isempty(b_nA_); path_gen2 = path_gen; path_gen2(b_nA_+2) = '0'; path_gen2(b_nA_+3) = 'x'; clear_covariates_flag=1; lowrank_vs_diffexp_flag= 0; end;
b_nAx = strfind(path_gen,'_nAx'); if ~isempty(b_nAx); path_gen2 = path_gen; path_gen2(b_nAx+2) = '0'; lowrank_vs_diffexp_flag= 0; end;
if exist(sprintf('%s%s.mat',path_base,gen_fname),'file'); 
load(sprintf('%s%s.mat',path_base,gen_fname),'npats','ngenes','Pnames','Gnames','g_ij','D','d_ij','x_ij','cov_mat','cov_cat'); 
elseif ~exist(sprintf('%s%s.mat',path_base,gen_fname),'file') & ~isempty(path_gen2) & exist(path_gen2,'file'); 
if disp_flag; disp(sprintf(' %% loading %s rather than %s',path_gen2,path_gen)); end;
load(path_gen2,'npats','ngenes','Pnames','Gnames','g_ij','D','d_ij','x_ij','cov_mat','cov_cat');
else disp(sprintf(' %% Warning! %s does not exist!',path_gen));
end;% if ~exist;
ngenes_sub = length(g_ij); if (lowrank_vs_diffexp_flag==0); MATRIX_ORIG = D; end;
if clear_covariates_flag; if disp_flag; disp(sprintf(' %% reducing ncov from %d to 0',size(cov_mat,2))); end; cov_mat = zeros(npats,0); cov_cat = ones(npats,1); end;
if load_B_flag; load(path_gen2,'B'); if disp_flag; disp(sprintf(' %% loading matrix B [%d-x-%d]',size(B))); end; if size(B,2)==npats; B = transpose(B); end; end;
if disp_flag; disp(sprintf(' %% clear_covariates_flag %d load_B_flag %d lowrank_vs_diffexp_flag %d',clear_covariates_flag,load_B_flag,lowrank_vs_diffexp_flag)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
