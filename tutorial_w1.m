function tutorial_w1(path_base,gen_fname,gen_flag,prm_flag,n_to_find,frwd_vs_back,nthreads);
% function tutorial_w1();
% designed to read and process data such as GSE17895_n0x.mat ;
% try with: 
%{
  path_base = '/home/myname/dir_bicluster/dir_GSE17536/'; gen_fname = 'GSE17536_n1x'; 
  gen_flag = 1; prm_flag = 0; n_to_find = 1; frwd_vs_back = 'frwd'; nthreads = 2;  
  tutorial_w1(path_base,gen_fname,gen_flag,prm_flag,n_to_find,frwd_vs_back,nthreads);

  path_base = '/data/rangan/dir_bcc/dir_tutorial_081915/dir_GSE17536/'; gen_fname = 'GSE17536_n0x'; 
  gen_flag = 1; prm_flag = 0; n_to_find = 6; frwd_vs_back = 'frwd'; nthreads = 2;  
  tutorial_w1(path_base,gen_fname,gen_flag,prm_flag,n_to_find,frwd_vs_back,nthreads);

  %}

bitj=16; % this parameter is used in the compression routine 'tutorial_binary_compress.m' 
check_aggressive = 0;
prefix_base = sprintf('%s',gen_fname); prefix_use=prefix_base; 
path_pre_prm = sprintf('%sdir_%s_prm/',path_base,prefix_base);
if (prm_flag > 0 & ~exist(path_pre_prm,'dir')); command_string = sprintf('mkdir %s;',path_pre_prm); system(command_string); end;

if gen_flag | ~exist(sprintf('%s%s_gen.mat',path_base,prefix_base),'file');

tutorial_w1_setup;
Tall_n = [ones(npats,1) , cov_mat]; 
Vall_n = [zeros(1,ngenes)]; Vall_n(g_ij)=1; Vall_t = transpose(Vall_n);
if (lowrank_vs_diffexp_flag); tutorial_binary_compress(bitj,D>0,sprintf('%s%s_Aall_n.b16',path_base,prefix_base)); tutorial_binary_compress(bitj,transpose(D>0),sprintf('%s%s_Aall_t.b16',path_base,prefix_base));
 else save(sprintf('%s%s_Aall_n.mat',path_base,prefix_base),'MATRIX_ORIG'); end;
tutorial_binary_compress(bitj,Tall_n>0,sprintf('%s%s_Tall_n.b16',path_base,prefix_base));
tutorial_binary_compress(bitj,transpose(Tall_n)>0,sprintf('%s%s_Tall_t.b16',path_base,prefix_base));
tutorial_binary_compress(bitj,Vall_n>0,sprintf('%s%s_Vall_n.b16',path_base,prefix_base));
tutorial_binary_compress(bitj,transpose(Vall_n)>0,sprintf('%s%s_Vall_t.b16',path_base,prefix_base));
tutorial_binary_compress(bitj,ones(1,1),sprintf('%s%s_Vall_n_rind.b16',path_base,prefix_base));
tutorial_binary_compress(bitj,Vall_t,sprintf('%s%s_Ainc_n_cind.b16',path_base,prefix_base));
save(sprintf('%s%s_gen.mat',path_base,prefix_base),'bitj','gen_fname','npats','ngenes','ngenes_sub','g_ij','d_ij','x_ij','cov_mat','cov_cat','load_B_flag','clear_covariates_flag','lowrank_vs_diffexp_flag');

end;%if gen_flag;

if (n_to_find);

load_flag=1;
if load_flag;
load(sprintf('%s%s_gen.mat',path_base,prefix_base),'bitj','gen_fname','npats','ngenes','ngenes_sub','g_ij','d_ij','x_ij','cov_mat','cov_cat','load_B_flag','clear_covariates_flag','lowrank_vs_diffexp_flag');
end;%if load_flag;

prm_pre = 'prm'; if (prm_flag>0); prefix_use = sprintf('%s_%s%.3d',prefix_base,prm_pre,prm_flag); end;
for ns=max(1,min(cov_cat)):max(cov_cat);
std_d_ij{ns} = intersect(d_ij,find(cov_cat==ns));
std_x_ij{ns} = intersect(x_ij,find(cov_cat==ns));
end;%for ns=max(1,min(cov_cat)):max(cov_cat);

prm=[];
if (prm_flag);
% first we shuffle the random number generator, in case matlab starts with the same seed each time ;
rng(prm_flag); 
% now we shuffle labels within covariate classes;
for ns=max(1,min(cov_cat)):max(cov_cat);
disp(sprintf(' %% permuting labels within cov_cat %d',ns));
std_d_ij{ns} = intersect(d_ij,find(cov_cat==ns));
std_x_ij{ns} = intersect(x_ij,find(cov_cat==ns));
if (length(std_d_ij{ns})>0 & length(std_x_ij{ns})>0);
std_z_ij = [std_d_ij{ns} ; std_x_ij{ns}];
prm = randperm(length(std_z_ij));
std_d_ij{ns} = std_z_ij(prm(1:length(std_d_ij{ns})));
std_x_ij{ns} = std_z_ij(prm(length(std_d_ij{ns}) + (1:length(std_x_ij{ns}))));
end;%if (length(std_d_ij{ns})>0 & length(std_x_ij{ns})>0);
end;%for ns=max(1,min(cov_cat)):max(cov_cat);
end;% if (prm_flag);

% renormalize data if appropriate ;
renormalize_flag=0;
if prm_flag & load_B_flag ;
tutorial_w1_setup;
renormalize_flag=1;
E_tmp = B; D_tmp = zeros(size(E_tmp));
for ng = 1:ngenes;
if (mod(ng,1024)==0); disp(sprintf(' %% renormalizing gene %d out of %d',ng,ngenes)); end;
for ns=max(1,min(cov_cat)):max(cov_cat);
tmp = E_tmp(std_d_ij{ns},ng); tmp = (tmp-median(tmp))/std(tmp); D_tmp(std_d_ij{ns},ng) = tmp;
tmp = E_tmp(std_x_ij{ns},ng); tmp = (tmp-median(tmp))/std(tmp); D_tmp(std_x_ij{ns},ng) = tmp;
end;%for ns=max(1,min(cov_cat)):max(cov_cat);
end;%for ng = 1:ngenes;
path_use = path_pre_prm; 
if (lowrank_vs_diffexp_flag); tutorial_binary_compress(bitj,D_tmp>0,sprintf('%s%s_Aall_n.b16',path_use,prefix_use)); tutorial_binary_compress(bitj,transpose(D_tmp>0),sprintf('%s%s_Aall_t.b16',path_use,prefix_use));
 else MATRIX_ORIG = D_tmp; save(sprintf('%s%s_Aall_n.mat',path_use,prefix_use),'MATRIX_ORIG'); end;
end;% if permutation; 

disp(' %% Now we set up the case and control indices.');
Ainc_n_rij=[]; Zinc_n_rij=[];
for ns=max(1,min(cov_cat)):max(cov_cat);
Ainc_n_rij = union(Ainc_n_rij,std_d_ij{ns});
Zinc_n_rij = union(Zinc_n_rij,std_x_ij{ns});
end;%for ns=max(1,min(cov_cat)):max(cov_cat);

disp(' %% Now we set up the covariate indices.');
Tinc_n = [ones(length(Ainc_n_rij),1) , cov_mat(Ainc_n_rij,:)]; 
Sinc_n = [ones(length(Zinc_n_rij),1) , cov_mat(Zinc_n_rij,:)];

if (prm_flag>0); path_use = path_pre_prm; else path_use = path_base; end;

Ainc_n_rind = zeros(npats,1); Ainc_n_rind(Ainc_n_rij)=1;
tutorial_binary_compress(bitj,Ainc_n_rind,sprintf('%s%s_%s_Ainc_n_rind.b16',path_use,prefix_use,frwd_vs_back));
Zinc_n_rind = zeros(npats,1); Zinc_n_rind(Zinc_n_rij)=1;
tutorial_binary_compress(bitj,Zinc_n_rind,sprintf('%s%s_%s_Zinc_n_rind.b16',path_use,prefix_use,frwd_vs_back));
Tinc_n_cind = ones(1+size(cov_mat,2),1); 
tutorial_binary_compress(bitj,Tinc_n_cind,sprintf('%s%s_Tinc_n_cind.b16',path_use,prefix_use));

Ainc_n_rows = length(Ainc_n_rij);
Ainc_n_cols = ngenes;
Zinc_n_rows = length(Zinc_n_rij);
Zinc_n_cols = ngenes;
Tinc_n_cols = 1+size(cov_mat,2);

save(sprintf('%s%s_ij.mat',path_use,prefix_use),'bitj','gen_fname','ngenes','ngenes_sub','g_ij','npats','std_d_ij','std_x_ij','Ainc_n_rij','Zinc_n_rij','Tinc_n','Sinc_n','Ainc_n_rind','Zinc_n_rind','Tinc_n_cind','Ainc_n_rows','Ainc_n_cols','Zinc_n_rows','Zinc_n_cols','Tinc_n_cols');

if (prm_flag>0); % detecting dominant bicluster only ;

disp(sprintf(' %% Now we create a %s-input file for use with the code.',frwd_vs_back));
if (strcmp(frwd_vs_back,'frwd')); Dchar = 'A'; Xchar = 'Z'; elseif (strcmp(frwd_vs_back,'back')); Dchar = 'Z'; Xchar = 'A'; end;
fname_in = sprintf('%s%s_%s.in',path_use,prefix_use,frwd_vs_back);
fp=fopen(fname_in,'w');
fprintf(fp,'GLOBAL_verbose= 0;\n');
fprintf(fp,'GLOBAL_thread_count= %d;\n',nthreads); % choose number of threads;
fprintf(fp,'GLOBAL_CFILTER_ERRCHECK= 0;\n'); % set to 1 to doublecheck for errors midrun (slows down code immensely) ;
fprintf(fp,'GLOBAL_CFILTER_SCOREOUT= 0;\n'); % set to 1 to output scores at each iteration (generates massive files) ;
fprintf(fp,'GLOBAL_CFILTER_AGGRESSIVE= %0.2f;\n',check_aggressive); % set to 0 to eliminate rows/cols one by one ;
fprintf(fp,'GLOBAL_force_kr= 0;\n');
fprintf(fp,'GLOBAL_force_wk= 0;\n');
fprintf(fp,'GLOBAL_force_xv= +1;\n');
if renormalize_flag==0;
if (lowrank_vs_diffexp_flag); fprintf(fp,'A_n_name= %s%s_Aall_n.b16;\n',path_base,prefix_base); fprintf(fp,'A_t_name= %s%s_Aall_t.b16;\n',path_base,prefix_base);
 else fprintf(fp,'A_n_name= %s%s_Aall_n.mat;\n',path_base,prefix_base); end;
elseif renormalize_flag==1;
if (lowrank_vs_diffexp_flag); fprintf(fp,'A_n_name= %s%s_Aall_n.b16;\n',path_use,prefix_use); fprintf(fp,'A_t_name= %s%s_Aall_t.b16;\n',path_use,prefix_use);
 else fprintf(fp,'A_n_name= %s%s_Aall_n.mat;\n',path_use,prefix_use); end;
end;% if renormalize_flag;
fprintf(fp,'A_n_rows= %d;\n',npats);
fprintf(fp,'A_n_cols= %d;\n',ngenes);
fprintf(fp,'A_n_rind= %s%s_%s_%cinc_n_rind.b16;\n',path_use,prefix_use,frwd_vs_back,Dchar);
fprintf(fp,'A_n_cind= %s%s_Ainc_n_cind.b16;\n',path_base,prefix_base);
fprintf(fp,'Z_n_name= %s%s_Aall_n.b16;\n',path_base,prefix_base);
fprintf(fp,'Z_t_name= %s%s_Aall_t.b16;\n',path_base,prefix_base);
fprintf(fp,'Z_n_rows= %d;\n',npats);
fprintf(fp,'Z_n_rind= %s%s_%s_%cinc_n_rind.b16;\n',path_use,prefix_use,frwd_vs_back,Xchar);
fprintf(fp,'T_n_name= %s%s_Tall_n.b16;\n',path_base,prefix_base);
fprintf(fp,'T_t_name= %s%s_Tall_t.b16;\n',path_base,prefix_base);
fprintf(fp,'T_n_cols= %d;\n',Tinc_n_cols);
fprintf(fp,'T_n_cind= %s%s_Tinc_n_cind.b16;\n',path_use,prefix_use);
fprintf(fp,'S_n_name= %s%s_Tall_n.b16;\n',path_base,prefix_base);
fprintf(fp,'S_t_name= %s%s_Tall_t.b16;\n',path_base,prefix_base);
fprintf(fp,'V_n_name= %s%s_Vall_n.b16;\n',path_base,prefix_base);
fprintf(fp,'V_t_name= %s%s_Vall_t.b16;\n',path_base,prefix_base);
fprintf(fp,'V_n_rows= %d;\n',1);
fprintf(fp,'V_n_rind= %s%s_Vall_n_rind.b16;\n',path_base,prefix_base);
fprintf(fp,'GLOBAL_out_name= %s%s_%s_out;\n',path_use,prefix_use,frwd_vs_back);
fprintf(fp,'GLOBAL_prm_flag= %d;\n',prm_flag);
fprintf(fp,'END= 0;\n');
fclose(fp);
disp(' %% We run the matlab-code...');
if (lowrank_vs_diffexp_flag==1); tutorial_lakcluster_1(sprintf('%s',fname_in)); elseif (lowrank_vs_diffexp_flag==2); tutorial_lakcluster_2(sprintf('%s',fname_in)); else if clear_covariates_flag; tutorial_dexcluster_1(sprintf('%s',fname_in)); elseif (~clear_covariates_flag); tutorial_dexcluster_2(sprintf('%s',fname_in)); end; end;
if renormalize_flag==1;
if (lowrank_vs_diffexp_flag); delete(sprintf('%s%s_Aall_n.b16',path_use,prefix_use)); delete(sprintf('%s%s_Aall_t.b16',path_use,prefix_use));
 else delete(sprintf('%s%s_Aall_n.mat',path_use,prefix_use)); end;
end;% if renormalize_flag;
delete(sprintf('%s%s_%s.in',path_use,prefix_use,frwd_vs_back));
delete(sprintf('%s%s_%s_Ainc_n_rind.b16',path_use,prefix_use,frwd_vs_back));
delete(sprintf('%s%s_%s_Zinc_n_rind.b16',path_use,prefix_use,frwd_vs_back));
delete(sprintf('%s%s_Tinc_n_cind.b16',path_use,prefix_use));
delete(sprintf('%s%s_ij.mat',path_use,prefix_use));
delete(sprintf('%s%s_%s_out_xdrop.txt',path_use,prefix_use,frwd_vs_back));

elseif (prm_flag==0); % detecting biclusters, removing, and rerunning ;

nrun=0;
while (nrun<=n_to_find-1);
if (nrun==0); 
path_plus_prefix_use = sprintf('%s%s_%s',path_use,prefix_use,frwd_vs_back); 
path_plus_txt_use = sprintf('%sdir_txt/',path_use); if (~exist(path_plus_txt_use,'dir')); mkdir(path_plus,'dir_txt'); end;
path_plus_txt_plus_prefix_use = sprintf('%s%s_%s',path_plus_txt_use,prefix_use,frwd_vs_back); 
 else path_plus_prefix_use = sprintf('%s%s_%s_%d',path_use,prefix_use,frwd_vs_back,nrun); end;%if (nrun==0);

if (strcmp(frwd_vs_back,'frwd')); Dchar = 'A'; Xchar = 'Z'; elseif (strcmp(frwd_vs_back,'back')); Dchar = 'Z'; Xchar = 'A'; end;
fname_in = sprintf('%s.in',path_plus_prefix_use);
disp(sprintf(' %% Now we create a %s-input file %s (nrun %d) for use with the code.',frwd_vs_back,fname_in,nrun));
fp=fopen(fname_in,'w');
fprintf(fp,'GLOBAL_verbose= 0;\n');
fprintf(fp,'GLOBAL_thread_count= %d;\n',nthreads); % choose number of threads;
fprintf(fp,'GLOBAL_CFILTER_ERRCHECK= 0;\n'); % set to 1 to doublecheck for errors midrun (slows down code immensely) ;
fprintf(fp,'GLOBAL_CFILTER_SCOREOUT= 0;\n'); % set to 1 to output scores at each iteration (generates massive files) ;
fprintf(fp,'GLOBAL_CFILTER_AGGRESSIVE= %0.2f;\n',check_aggressive); % set to 0 to eliminate rows/cols one by one ;
fprintf(fp,'GLOBAL_force_kr= 0;\n');
fprintf(fp,'GLOBAL_force_wk= 0;\n');
fprintf(fp,'GLOBAL_force_xv= +1;\n');
if (lowrank_vs_diffexp_flag); fprintf(fp,'A_n_name= %s%s_Aall_n.b16;\n',path_base,prefix_base); fprintf(fp,'A_t_name= %s%s_Aall_t.b16;\n',path_base,prefix_base);
 else fprintf(fp,'A_n_name= %s%s_Aall_n.mat;\n',path_base,prefix_base); end;
fprintf(fp,'A_n_rows= %d;\n',npats);
fprintf(fp,'A_n_cols= %d;\n',ngenes);
fprintf(fp,'A_n_rind= %s%s_%s_%cinc_n_rind.b16;\n',path_use,prefix_use,frwd_vs_back,Dchar);
fprintf(fp,'A_n_cind= %s%s_Ainc_n_cind.b16;\n',path_base,prefix_base);
fprintf(fp,'Z_n_name= %s%s_Aall_n.b16;\n',path_base,prefix_base);
fprintf(fp,'Z_t_name= %s%s_Aall_t.b16;\n',path_base,prefix_base);
fprintf(fp,'Z_n_rows= %d;\n',npats);
fprintf(fp,'Z_n_rind= %s%s_%s_%cinc_n_rind.b16;\n',path_use,prefix_use,frwd_vs_back,Xchar);
fprintf(fp,'T_n_name= %s%s_Tall_n.b16;\n',path_base,prefix_base);
fprintf(fp,'T_t_name= %s%s_Tall_t.b16;\n',path_base,prefix_base);
fprintf(fp,'T_n_cols= %d;\n',Tinc_n_cols);
fprintf(fp,'T_n_cind= %s%s_Tinc_n_cind.b16;\n',path_use,prefix_use);
fprintf(fp,'S_n_name= %s%s_Tall_n.b16;\n',path_base,prefix_base);
fprintf(fp,'S_t_name= %s%s_Tall_t.b16;\n',path_base,prefix_base);
fprintf(fp,'V_n_name= %s%s_Vall_n.b16;\n',path_base,prefix_base);
fprintf(fp,'V_t_name= %s%s_Vall_t.b16;\n',path_base,prefix_base);
fprintf(fp,'V_n_rows= %d;\n',1);
fprintf(fp,'V_n_rind= %s%s_Vall_n_rind.b16;\n',path_base,prefix_base);
fprintf(fp,'GLOBAL_out_name= %s_out;\n',path_plus_prefix_use);
fprintf(fp,'GLOBAL_prm_flag= %d;\n',prm_flag);
if (nrun>0);
fprintf(fp,'A_n_repl_num= %d;\n',nrun);
fprintf(fp,'A_n_repl_name= ');
for nrun_index=0:nrun-1;
if (nrun_index==0); fprintf(fp,'%s%s_%s_bc.txt',path_use,prefix_use,frwd_vs_back); 
 else; fprintf(fp,'%s%s_%s_%d_bc.txt',path_use,prefix_use,frwd_vs_back,nrun_index); end;% if (nrun_index==0);
if (nrun_index<nrun-1); fprintf(fp,', '); else fprintf(fp,';\n'); end;
end;%for nrun_index=0:nrun-1;
end;%if (nrun>0);
fprintf(fp,'END= 0;\n');
fclose(fp);
disp(' %% We run the matlab-code...');
if (lowrank_vs_diffexp_flag==1); tutorial_lakcluster_1(sprintf('%s',fname_in)); elseif (lowrank_vs_diffexp_flag==2); tutorial_lakcluster_2(sprintf('%s',fname_in)); else if clear_covariates_flag; tutorial_dexcluster_1(sprintf('%s',fname_in)); elseif ~clear_covariates_flag; tutorial_dexcluster_2(sprintf('%s',fname_in)); end; end;

% calculating loopS ;
if (lowrank_vs_diffexp_flag==1);
if strcmp(frwd_vs_back,'frwd'); tutorial_lakcluster_1(sprintf('%s',fname_in),1,min(512,Ainc_n_rows),min(4096,length(g_ij))); end;
if strcmp(frwd_vs_back,'back'); tutorial_lakcluster_1(sprintf('%s',fname_in),1,min(512,Zinc_n_rows),min(4096,length(g_ij))); end;
 elseif (lowrank_vs_diffexp_flag==2);
if strcmp(frwd_vs_back,'frwd'); tutorial_lakcluster_2(sprintf('%s',fname_in),1,min(512,Ainc_n_rows),min(4096,length(g_ij))); end;
if strcmp(frwd_vs_back,'back'); tutorial_lakcluster_2(sprintf('%s',fname_in),1,min(512,Zinc_n_rows),min(4096,length(g_ij))); end;
elseif (lowrank_vs_diffexp_flag==0) & (clear_covariates_flag);
if strcmp(frwd_vs_back,'frwd'); tutorial_dexcluster_1(sprintf('%s',fname_in),1,min(512,Ainc_n_rows),min(4096,length(g_ij))); end;
if strcmp(frwd_vs_back,'back'); tutorial_dexcluster_1(sprintf('%s',fname_in),1,min(512,Zinc_n_rows),min(4096,length(g_ij))); end;
elseif (lowrank_vs_diffexp_flag==0) & (~clear_covariates_flag); 
if strcmp(frwd_vs_back,'frwd'); tutorial_dexcluster_2(sprintf('%s',fname_in),1,min(512,Ainc_n_rows),min(4096,length(g_ij))); end;
if strcmp(frwd_vs_back,'back'); tutorial_dexcluster_2(sprintf('%s',fname_in),1,min(512,Zinc_n_rows),min(4096,length(g_ij))); end;
end;% which algorithm to use for loopS ;
tutorial_w1_setup;
load(sprintf('%s%s_gen.mat',path_base,prefix_base),'bitj','gen_fname','npats','ngenes','ngenes_sub','g_ij','d_ij','x_ij','cov_mat','cov_cat','load_B_flag','clear_covariates_flag','lowrank_vs_diffexp_flag');
load(sprintf('%s%s_ij.mat',path_use,prefix_use),'bitj','gen_fname','ngenes','ngenes_sub','g_ij','npats','std_d_ij','std_x_ij','Ainc_n_rij','Zinc_n_rij','Tinc_n','Sinc_n','Ainc_n_rind','Zinc_n_rind','Tinc_n_cind','Ainc_n_rows','Ainc_n_cols','Zinc_n_rows','Zinc_n_cols','Tinc_n_cols');
load(sprintf('%s_out_loopS.mat',path_plus_prefix_use),'ncbins','nrbins','ncmax','nrmax','nc_ra','nr_ra','LR_ra','LC_ra','Lx_ra','LR2_ra','LC2_ra','Lx2_ra');

plot_flag=1;
if (plot_flag);
figure; cla;
subplot(1,5,[1,2]); imagesc(LR2_ra,[0,1]); colorbar; title('average row score'); set(gca,'Ytick',[],'Xtick',[]); xlabel('columns n'); ylabel('rows m');
subplot(1,5,[3,4]); imagesc(LC2_ra,[0,1]); colorbar; title('average col score'); set(gca,'Ytick',[],'Xtick',[]); xlabel('columns n'); ylabel('rows m');
subplot(1,5,5); stairs(Lx2_ra(end:-1:1,1),1:size(Lx2_ra,1)); xlim([0,max(Lx2_ra(:))+1]); ylim([1,size(Lx2_ra,1)]); set(gca,'Ytick',[]); title('# covariates remaining');
xlabel('covariate-categories remaining'); ylabel('rows m');
print('-depsc',sprintf('%s_loopS.eps',path_plus_prefix_use));
print('-djpeg',sprintf('%s_loopS.jpg',path_plus_prefix_use));
end;%if (plot_flag);
if (strcmp(frwd_vs_back,'frwd')); Dinc_n_rij = Ainc_n_rij; Xinc_n_rij = Zinc_n_rij; end;
if (strcmp(frwd_vs_back,'back')); Dinc_n_rij = Zinc_n_rij; Xinc_n_rij = Ainc_n_rij; end;
out_xdrop = textread(sprintf('%s_out_xdrop.txt',path_plus_prefix_use));
rdrop = out_xdrop(:,1); cdrop = out_xdrop(:,2);
rij = rdrop(find(rdrop>=0)); rij = [rij(:) ; setdiff(Dinc_n_rij(:)-1,rij(:))]; rij = rij(end:-1:1);
cij = cdrop(find(cdrop>=0)); cij = [cij(:) ; setdiff(g_ij(:)-1,cij(:))]; cij = cij(end:-1:1);
rthr = 0.7; ar_min = 0.1 ; ar_max = 1.0/ar_min ;
[output,ik_out,ij_out] = tutorial_loopS_threshold(LR2_ra,Lx2_ra,rthr,ar_min,ar_max);

if (output>0);
plot_flag=1;
ncij = 1:nc_ra(ik_out); nrij = 1:nr_ra(ij_out);
rij_sub = rij(nrij); cij_sub = cij(ncij);
pgap = 32;ggap = 4;
ncovs=size(cov_mat,2);cdup=zeros(ncovs,ncovs*ggap);for nl=0:ncovs-1;cdup(1+nl,1+nl*ggap+(1:ggap))=1;end;
cov_tmp = cov_mat*cdup;
cmap_autumn = colormap('autumn'); cmap_autumn(1,:) = [1,1,1]; cmap_autumn(end-1,:) = 0.85*[1,1,1]; cmap_autumn(end,:) = 0.00*[1,1,1]; clen = size(cmap_autumn,1); colormap(cmap_autumn);
cbot=-1.5; ctop=1.5; cvals = linspace(cbot,ctop,clen); cbotl=cvals(2); ctopl = cvals(end-2); cdiff = cvals(end)-cvals(end-1); 
Dtmp = D(1+rij_sub ,1+cij_sub) ; cov_tmp_D = cov_tmp(1+rij_sub ,:);
Xtmp = D(Zinc_n_rij,1+cij_sub) ; cov_tmp_X = cov_tmp(Zinc_n_rij,:);
subplot(1,2,1);
imagesc([...
  min(ctopl,max(cbotl,Dtmp)) , cbot*ones(length(1+rij_sub),ggap) , ((ctopl+cdiff+cdiff*cov_tmp_D)) ; ...
  cbot*ones(pgap,length(cij_sub)) , cbot*ones(pgap,ggap) , cbot*ones(pgap,size(cov_tmp,2)); ...
  min(ctopl,max(cbotl,Xtmp)) , cbot*ones(length( x_ij),ggap) , ((ctopl+cdiff+cdiff*cov_tmp_X)) ; ...
	 ],[cbot,ctop]);
set(gca,'Xtick',[],'Ytick',[]); 
title(sprintf('%s_bc.txt',path_plus_prefix_use),'interpreter','none');
subplot(1,2,2);
[prA,pc,prX] = tutorial_hcluster(Dtmp,Xtmp); Dtmp = Dtmp(prA,pc); Xtmp = Xtmp(prX,pc); cov_tmp_D = cov_tmp_D(prA,:); cov_tmp_X = cov_tmp_X(prX,:);
imagesc([...
  min(ctopl,max(cbotl,Dtmp)) , cbot*ones(length(1+rij_sub),ggap) , ((ctopl+cdiff+cdiff*cov_tmp_D)) ; ...
  cbot*ones(pgap,length(cij_sub)) , cbot*ones(pgap,ggap) , cbot*ones(pgap,size(cov_tmp,2)); ...
  min(ctopl,max(cbotl,Xtmp)) , cbot*ones(length( x_ij),ggap) , ((ctopl+cdiff+cdiff*cov_tmp_X)) ; ...
	 ],[cbot,ctop]);
set(gca,'Xtick',[],'Ytick',[]); 
title(sprintf('%s_bc.txt (rearranged)',path_plus_prefix_use),'interpreter','none');
if (plot_flag);
print('-depsc',sprintf('%s_bc.eps',path_plus_prefix_use));
print('-djpeg',sprintf('%s_bc.jpg',path_plus_prefix_use));
end;%if (plot_flag);
rij_sub = rij(nrij); cij_sub = cij(ncij);
bc_fname = sprintf('%s_bc.txt',path_plus_prefix_use);
fid = fopen(bc_fname,'w');
for nl=1:length(rij_sub); fprintf(fid,'%d %d\n',rij_sub(nl),-1); end;%for nl=1:length(rij_sub);
for nl=1:length(cij_sub); fprintf(fid,'%d %d\n',-1,cij_sub(nl)); end;%for nl=1:length(cij_sub);
fclose(fid);
rij_sub = rij(nrij);
plist_fname = sprintf('%s_plist.txt',path_plus_txt_plus_prefix_use); 
% Warning! plist indices formed by [iv,hdr_ij,dat_ij] = intersect(hdr_val,dat_val,'stable'); If hdr_val and dat_val are not organized correctly, these plist indices could be misleading. ;
fid = fopen(plist_fname,'w');
for nl=1:length(rij_sub);
tmp = sprintf('%d %d',Pnames(1+rij_sub(nl)),1+rij_sub(nl));
for nl2=1:size(cov_mat,2);
tmp = sprintf('%s %d',tmp,cov_mat(1+rij(nl),nl2));
end;%for nl2=1:size(cov_mat,2);
tmp = sprintf('%s %d',tmp,Lx2_ra(nl,end));
fprintf(fid,'%s\n',tmp);
end;%for nl=1:length(rij_sub);
fclose(fid);
cij_sub = cij(ncij);
g_id_sub = Gnames(1+cij_sub);
glist_fname = sprintf('%s_glist.txt',path_plus_txt_plus_prefix_use);
fid = fopen(glist_fname,'w');
disp(sprintf('found %d genes (across %d patients): see %s for list',max(ncij),max(nrij),glist_fname));
for nl=1:length(g_id_sub);
fprintf(fid,'%d\n',g_id_sub(nl));
end;%for nl=1:length(g_id_sub);
fclose(fid);
gname_fname = sprintf('%s_gname.txt',path_plus_txt_plus_prefix_use);
tutorial_gene_number_to_name(path_base,g_id_sub,gname_fname);
%genri_fname = sprintf('%s_genri.txt',path_plus_prefix_use);
%command_string = sprintf('java seek.GeneEnrichmentTest -l /data/rangan/dir_bcc/dir_code_012015/hierarchical/enrichment/ -m 2 %s -s original -t %d >! %s ;',glist_fname,length(g_id_sub),genri_fname);
%disp(command_string); system(command_string);
%genri_fname = sprintf('%s_gslim.txt',path_plus_prefix_use);
%command_string = sprintf('java seek.GeneEnrichmentTest -l /data/rangan/dir_bcc/dir_code_012015/hierarchical/enrichment/ -m 2 %s -s experimental_bp_slim -t %d >! %s ;',glist_fname,length(g_id_sub),genri_fname);
%disp(command_string); system(command_string);
nrun = nrun + 1;
 else;%if (output==0);
disp('helper_loopS_to_genes failed to find bicluster, stopping search for more biclusters. ');
nrun=nn_to_find;
end;%if (output>0);

end;%while (nrun<=n_to_find-1);

end;%if (prm_flag==0); % detecting biclusters, removing, and rerunning ;

end;%if (n_to_find);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

