% This generates an example and tests tutorial_lakcluster_1.m ;

disp(sprintf(' '));
disp(sprintf(' This function tests tutorial_lakcluster_1.'));
rng(1);
path_base = './dir_tutorial_lakcluster_test/'; if (~exist(path_base,'dir')); mkdir(path_base); end;
prefix_base = 'tutorial_lakcluster_1_test';
bitj = 16; prm_flag = 0; nthreads = 2; check_aggressive = 0.025;
npats = 1024; ngenes = 2048; g_ij = 1:ngenes;
Aall_n = randn(npats,ngenes);
ptmp = randperm(npats); d_ij = ptmp(1:floor(npats/2)); x_ij = setdiff(1:npats,d_ij);
ptmp = randperm(length(d_ij)); c1d__on = ptmp(1:floor(length(d_ij)/2)); c1d_off = setdiff(1:length(d_ij),c1d__on); c1d_ord = [c1d__on,c1d_off]; [tmp,c1d_ord_s] = sort(c1d_ord); [tmp,c1d_ord_r] = sort(c1d_ord_s);
ptmp = randperm(length(x_ij)); c1x__on = ptmp(1:floor(length(x_ij)/2)); c1x_off = setdiff(1:length(x_ij),c1x__on); c1x_ord = [c1x__on,c1x_off]; [tmp,c1x_ord_s] = sort(c1x_ord); [tmp,c1x_ord_r] = sort(c1x_ord_s);
ptmp = randperm(length(d_ij)); c2d__on = ptmp(1:floor(length(d_ij)/2)); c2d_off = setdiff(1:length(d_ij),c2d__on); c2d_ord = [c2d__on,c2d_off]; [tmp,c2d_ord_s] = sort(c2d_ord); [tmp,c2d_ord_r] = sort(c2d_ord_s);
ptmp = randperm(length(x_ij)); c2x__on = ptmp(1:floor(length(x_ij)/2)); c2x_off = setdiff(1:length(x_ij),c2x__on); c2x_ord = [c2x__on,c2x_off]; [tmp,c2x_ord_s] = sort(c2x_ord); [tmp,c2x_ord_r] = sort(c2x_ord_s);
C0rows = floor(length(d_ij).^(0.65)); C0cols = floor(length(g_ij).^(0.65)); C0 = randn(C0rows,1)*randn(1,C0cols) + 0.01*randn(C0rows,C0cols); [prA,pc] = tutorial_hcluster(C0); C0=C0(prA,pc);
C1rows = floor(length(d_ij).^(0.70)); C1cols = floor(length(g_ij).^(0.70)); C1 = randn(C1rows,1)*randn(1,C1cols) + 0.01*randn(C1rows,C1cols); [prA,pc] = tutorial_hcluster(C1); C1=C1(prA,pc);
C2rows = floor(length(d_ij).^(0.70)); C2cols = floor(length(g_ij).^(0.70)); C2 = randn(C2rows,1)*randn(1,C2cols) + 0.01*randn(C2rows,C2cols); [prA,pc] = tutorial_hcluster(C2); C2=C2(prA,pc);
C3rows = 2*floor(length(d_ij).^(0.75)); C3cols = floor(length(g_ij).^(0.75)); C3 = randn(C3rows,1)*randn(1,C3cols) + 0.01*randn(C3rows,C3cols); 
[prA,pc] = tutorial_hcluster(C3(1:C3rows/2,:)); C3(1:C3rows/2,:)=C3(prA,pc);
[prA,pc] = tutorial_hcluster(C3(C3rows/2+(1:C3rows/2),:)); C3(C3rows/2 + (1:C3rows/2),:)=C3(C3rows/2 + (prA),pc);
p_d_ij = randperm(length(d_ij));[p_d_ij_v,p_d_ij_s] = sort(p_d_ij); [p_d_ij_v,p_d_ij_r] = sort(p_d_ij_s);
q_d_ij = randperm(length(d_ij));[q_d_ij_v,q_d_ij_s] = sort(q_d_ij); [q_d_ij_v,q_d_ij_r] = sort(q_d_ij_s);
q_x_ij = randperm(length(x_ij));[q_x_ij_v,q_x_ij_s] = sort(q_x_ij); [q_x_ij_v,q_x_ij_r] = sort(q_x_ij_s);
p0_g_ij = randperm(length(g_ij));[p0_g_ij_v,p0_g_ij_s] = sort(p0_g_ij); [p0_g_ij_v,p0_g_ij_r] = sort(p0_g_ij_s);
p1_g_ij = randperm(length(g_ij));[p1_g_ij_v,p1_g_ij_s] = sort(p1_g_ij); [p1_g_ij_v,p1_g_ij_r] = sort(p1_g_ij_s);
p2_g_ij = randperm(length(g_ij));[p2_g_ij_v,p2_g_ij_s] = sort(p2_g_ij); [p2_g_ij_v,p2_g_ij_r] = sort(p2_g_ij_s);
p3_g_ij = p0_g_ij(end:-1:1); p3_g_ij_r = p0_g_ij_r(end:-1:1);
Aall_n(d_ij(q_d_ij(1:C3rows/2)),g_ij(p3_g_ij(1:C3cols))) = C3(1:C3rows/2,:);
Aall_n(x_ij(q_x_ij(1:C3rows/2)),g_ij(p3_g_ij(1:C3cols))) = C3(C3rows/2 + (1:C3rows/2),:);
Aall_n(d_ij(c2d_ord(1:C2rows)),g_ij(p2_g_ij(1:C2cols))) = C2;
Aall_n(d_ij(c1d_ord(1:C1rows)),g_ij(p1_g_ij(1:C1cols))) = C1;
Aall_n(d_ij(p_d_ij(1:C0rows)),g_ij(p0_g_ij(1:C0cols))) = C0;
Aall_n = 2*(Aall_n>0)-1; % binarization step ;

disp(sprintf(' '));
disp(sprintf(' We start out by generating some data, comprising 1024 patients (512 cases and 512 controls) and 2048 genes. '));
disp(sprintf(' This data is also equipped with two categorical covariates (covariate-1 and covariate-2). '));
disp(sprintf(' Within this data we embed 4 biclusters: '));
disp(sprintf(' bc0 : this bicluster is case-specific, and is balanced across the covariates. This is the one we want to find. '));
disp(sprintf(' bc1 : this bicluster is also case-specific, but is biased with respect to covariate-1. We want to ignore this. '));
disp(sprintf(' bc2 : this bicluster is also case-specific, but is biased with respect to covariate-2. We want to ignore this. '));
disp(sprintf(' bc3 : this bicluster extends across both the cases and controls. We want to ignore this also. '));
disp(sprintf(' '));
disp(sprintf(' Now we plot the data four times, '));
disp(sprintf(' each time choosing the row/column ordering which reveals one of the embedded biclusters: '));
disp(sprintf(' subplots from left to right show bc0 --> bc3, with the case-matrix on top and the control-matrix below. '));
disp(sprintf(' Within each subplot the covariates are plotted as the broad gray/black columns on the right-hand side. '));
disp(sprintf(' Note that bc0 is the only bicluster which is both case-specific and balanced with respect to the covariates. '));
ncovs = 2;
cov_mat = zeros(npats,2); 
cov_mat(d_ij(c1d__on),1) = 1; cov_mat(x_ij(c1x__on),1) = 1;
cov_mat(d_ij(c2d__on),2) = 1; cov_mat(x_ij(c2x__on),2) = 1;
cov_cat = [ones(npats,1) , cov_mat]*transpose(2.^(0:ncovs));  
pgap = 32;ggap = 128;
cov_tmp = zeros(size(cov_mat,1),size(cov_mat,2)*ggap); for nc=1:size(cov_mat,2); cov_tmp(:,(nc-1)*ggap+(1:ggap)) = repmat(cov_mat(:,nc),1,ggap); end;%for nc=1:size(cov_mat,2);
cmap_autumn = colormap('autumn'); cmap_autumn(1,:) = [1,1,1]; cmap_autumn(end-1,:) = 0.85*[1,1,1]; cmap_autumn(end,:) = 0.00*[1,1,1]; clen = size(cmap_autumn,1); colormap(cmap_autumn);
cbot=-1.5; ctop=1.5; cvals = linspace(cbot,ctop,clen); cbotl=cvals(2); ctopl = cvals(end-2); cdiff = cvals(end)-cvals(end-1); 
subplot(1,4,1); 
Dtmp = Aall_n(d_ij(p_d_ij_r),g_ij(p0_g_ij_r)) ; cov_tmp_D = cov_tmp(d_ij(p_d_ij_r),:);
Xtmp = Aall_n(x_ij(:),g_ij(p0_g_ij_r)) ; cov_tmp_X = cov_tmp(x_ij(:),:);
colormap(cmap_autumn);
imagesc([...
	   min(ctopl,max(cbotl,Dtmp)) , cbot*ones(size(Dtmp,1),ggap) , ((ctopl+cdiff+cdiff*cov_tmp_D)) ; ...
  cbot*ones(pgap,ngenes) , cbot*ones(pgap,ggap) , cbot*ones(pgap,size(cov_tmp,2)); ...
	      min(ctopl,max(cbotl,Xtmp)) , cbot*ones(size(Xtmp,1),ggap) , ((ctopl+cdiff+cdiff*cov_tmp_X)) ; ...
	 ],[cbot,ctop]);
  set(gca,'Xtick',[],'Ytick',[]); title('case-only bicluster 0');
subplot(1,4,2); 
Dtmp = Aall_n(d_ij(c1d_ord_r),g_ij(p1_g_ij_r)) ; cov_tmp_D = cov_tmp(d_ij(c1d_ord_r),:);
Xtmp = Aall_n(x_ij(c1x_ord_r),g_ij(p1_g_ij_r)) ; cov_tmp_X = cov_tmp(x_ij(c1x_ord_r),:);
colormap(cmap_autumn);
imagesc([...
	   min(ctopl,max(cbotl,Dtmp)) , cbot*ones(size(Dtmp,1),ggap) , ((ctopl+cdiff+cdiff*cov_tmp_D)) ; ...
  cbot*ones(pgap,ngenes) , cbot*ones(pgap,ggap) , cbot*ones(pgap,size(cov_tmp,2)); ...
	      min(ctopl,max(cbotl,Xtmp)) , cbot*ones(size(Xtmp,1),ggap) , ((ctopl+cdiff+cdiff*cov_tmp_X)) ; ...
	 ],[cbot,ctop]);
  set(gca,'Xtick',[],'Ytick',[]); title('covariate-imbalanced bicluster 1');
subplot(1,4,3); 
Dtmp = Aall_n(d_ij(c2d_ord_r),g_ij(p2_g_ij_r)) ; cov_tmp_D = cov_tmp(d_ij(c2d_ord_r),:);
Xtmp = Aall_n(x_ij(c2x_ord_r),g_ij(p2_g_ij_r)) ; cov_tmp_X = cov_tmp(x_ij(c2x_ord_r),:);
colormap(cmap_autumn);
imagesc([...
	   min(ctopl,max(cbotl,Dtmp)) , cbot*ones(size(Dtmp,1),ggap) , ((ctopl+cdiff+cdiff*cov_tmp_D)) ; ...
  cbot*ones(pgap,ngenes) , cbot*ones(pgap,ggap) , cbot*ones(pgap,size(cov_tmp,2)); ...
	      min(ctopl,max(cbotl,Xtmp)) , cbot*ones(size(Xtmp,1),ggap) , ((ctopl+cdiff+cdiff*cov_tmp_X)) ; ...
	 ],[cbot,ctop]);
  set(gca,'Xtick',[],'Ytick',[]); title('covariate-imbalanced bicluster 2');
subplot(1,4,4); 
Dtmp = Aall_n(d_ij(q_d_ij_r),g_ij(p3_g_ij_r)) ; cov_tmp_D = cov_tmp(d_ij(q_d_ij_r),:);
Xtmp = Aall_n(x_ij(q_x_ij_r),g_ij(p3_g_ij_r)) ; cov_tmp_X = cov_tmp(x_ij(q_x_ij_r),:);
colormap(cmap_autumn);
imagesc([...
	   min(ctopl,max(cbotl,Dtmp)) , cbot*ones(size(Dtmp,1),ggap) , ((ctopl+cdiff+cdiff*cov_tmp_D)) ; ...
  cbot*ones(pgap,ngenes) , cbot*ones(pgap,ggap) , cbot*ones(pgap,size(cov_tmp,2)); ...
	      min(ctopl,max(cbotl,Xtmp)) , cbot*ones(size(Xtmp,1),ggap) , ((ctopl+cdiff+cdiff*cov_tmp_X)) ; ...
	 ],[cbot,ctop]);
  set(gca,'Xtick',[],'Ytick',[]); title('case-control bicluster 3');

disp(sprintf(' '));
disp(sprintf(' Paused (push a key). '));
pause();
disp(sprintf(' '));
disp(sprintf(' Now we set up the input file for the tutorial_lakcluster_1.m code. '));

Tall_n = [ones(npats,1) , cov_mat]; 
Vall_n = [zeros(1,ngenes)]; Vall_n(g_ij)=1; Vall_t = transpose(Vall_n);
tutorial_binary_compress(bitj,Aall_n>0,sprintf('%s%s_Aall_n.b16',path_base,prefix_base));
tutorial_binary_compress(bitj,transpose(Aall_n)>0,sprintf('%s%s_Aall_t.b16',path_base,prefix_base));
tutorial_binary_compress(bitj,Tall_n>0,sprintf('%s%s_Tall_n.b16',path_base,prefix_base));
tutorial_binary_compress(bitj,transpose(Tall_n)>0,sprintf('%s%s_Tall_t.b16',path_base,prefix_base));
tutorial_binary_compress(bitj,Vall_n>0,sprintf('%s%s_Vall_n.b16',path_base,prefix_base));
tutorial_binary_compress(bitj,transpose(Vall_n)>0,sprintf('%s%s_Vall_t.b16',path_base,prefix_base));
tutorial_binary_compress(bitj,ones(1,1),sprintf('%s%s_Vall_n_rind.b16',path_base,prefix_base));
tutorial_binary_compress(bitj,Vall_t,sprintf('%s%s_Ainc_n_cind.b16',path_base,prefix_base));
disp(' First we set up the case and control indices.');
for ns=max(1,min(cov_cat)):max(cov_cat);
std_d_ij{ns} = intersect(d_ij,find(cov_cat==ns));
std_x_ij{ns} = intersect(x_ij,find(cov_cat==ns));
end;%for ns=max(1,min(cov_cat)):max(cov_cat);
Ainc_n_rij=[]; Zinc_n_rij=[];
for ns=max(1,min(cov_cat)):max(cov_cat);
Ainc_n_rij = union(Ainc_n_rij,std_d_ij{ns});
Zinc_n_rij = union(Zinc_n_rij,std_x_ij{ns});
end;%for ns=max(1,min(cov_cat)):max(cov_cat);
disp(' Then we set up the covariate indices.');
Tinc_n = [ones(length(Ainc_n_rij),1) , cov_mat(Ainc_n_rij,:)]; 
Sinc_n = [ones(length(Zinc_n_rij),1) , cov_mat(Zinc_n_rij,:)];
path_use = path_base; prefix = prefix_base;
Ainc_n_rind = zeros(npats,1); Ainc_n_rind(Ainc_n_rij)=1;
tutorial_binary_compress(bitj,Ainc_n_rind,sprintf('%s%s_Ainc_n_rind.b16',path_use,prefix));
Zinc_n_rind = zeros(npats,1); Zinc_n_rind(Zinc_n_rij)=1;
tutorial_binary_compress(bitj,Zinc_n_rind,sprintf('%s%s_Zinc_n_rind.b16',path_use,prefix));
Tinc_n_cind = ones(1+size(cov_mat,2),1); 
tutorial_binary_compress(bitj,Tinc_n_cind,sprintf('%s%s_Tinc_n_cind.b16',path_use,prefix));
Ainc_n_rows = length(Ainc_n_rij);
Ainc_n_cols = ngenes;
Zinc_n_rows = length(Zinc_n_rij);
Zinc_n_cols = ngenes;
Tinc_n_cols = 1+size(cov_mat,2);
frwd_vs_back = 'frwd';
disp(sprintf(' Then we create a %s-input file for use with the code.',frwd_vs_back));
if (strcmp(frwd_vs_back,'frwd')); Dchar = 'A'; Xchar = 'Z'; elseif (strcmp(frwd_vs_back,'back')); Dchar = 'Z'; Xchar = 'A'; end;
fname_in = sprintf('%s%s_%s.in',path_use,prefix,frwd_vs_back);
fp=fopen(fname_in,'w');
fprintf(fp,'GLOBAL_verbose= 0;\n');
fprintf(fp,'GLOBAL_thread_count= %d;\n',nthreads); % choose number of threads;
fprintf(fp,'GLOBAL_CFILTER_ERRCHECK= 0;\n'); % set to 1 to doublecheck for errors midrun (slows down code immensely) ;
fprintf(fp,'GLOBAL_CFILTER_SCOREOUT= 0;\n'); % set to 1 to output scores at each iteration (generates massive files) ;
fprintf(fp,'GLOBAL_CFILTER_AGGRESSIVE= %0.2f;\n',check_aggressive); % set to 0 to eliminate rows/cols one by one ;
fprintf(fp,'GLOBAL_force_kr= 0;\n');
fprintf(fp,'GLOBAL_force_wk= 0;\n');
fprintf(fp,'GLOBAL_force_xv= +1;\n');
fprintf(fp,'A_n_name= %s%s_Aall_n.b16;\n',path_base,prefix_base);
fprintf(fp,'A_t_name= %s%s_Aall_t.b16;\n',path_base,prefix_base);
fprintf(fp,'A_n_rows= %d;\n',npats);
fprintf(fp,'A_n_cols= %d;\n',ngenes);
fprintf(fp,'A_n_rind= %s%s_%cinc_n_rind.b16;\n',path_use,prefix,Dchar);
fprintf(fp,'A_n_cind= %s%s_Ainc_n_cind.b16;\n',path_base,prefix_base);
fprintf(fp,'Z_n_name= %s%s_Aall_n.b16;\n',path_base,prefix_base);
fprintf(fp,'Z_t_name= %s%s_Aall_t.b16;\n',path_base,prefix_base);
fprintf(fp,'Z_n_rows= %d;\n',npats);
fprintf(fp,'Z_n_rind= %s%s_%cinc_n_rind.b16;\n',path_use,prefix,Xchar);
fprintf(fp,'T_n_name= %s%s_Tall_n.b16;\n',path_base,prefix_base);
fprintf(fp,'T_t_name= %s%s_Tall_t.b16;\n',path_base,prefix_base);
fprintf(fp,'T_n_cols= %d;\n',Tinc_n_cols);
fprintf(fp,'T_n_cind= %s%s_Tinc_n_cind.b16;\n',path_use,prefix);
fprintf(fp,'S_n_name= %s%s_Tall_n.b16;\n',path_base,prefix_base);
fprintf(fp,'S_t_name= %s%s_Tall_t.b16;\n',path_base,prefix_base);
fprintf(fp,'V_n_name= %s%s_Vall_n.b16;\n',path_base,prefix_base);
fprintf(fp,'V_t_name= %s%s_Vall_t.b16;\n',path_base,prefix_base);
fprintf(fp,'V_n_rows= %d;\n',1);
fprintf(fp,'V_n_rind= %s%s_Vall_n_rind.b16;\n',path_base,prefix_base);
fprintf(fp,'GLOBAL_out_name= %s%s_%s_out;\n',path_use,prefix,frwd_vs_back);
fprintf(fp,'GLOBAL_prm_flag= %d;\n',prm_flag);
fprintf(fp,'END= 0;\n');
fclose(fp);
disp(' Finally, we run the matlab-code...');
disp(sprintf(' '));
disp(sprintf(' This code will run through several iterations, removing rows and columns that are not '));
disp(sprintf(' obviously part of any bicluster, and retaining those that are. Hopefully, at the end, '));
disp(sprintf(' the rows and columns which are retained will be those corresponding to bicluster 0. '));
tutorial_lakcluster_1(sprintf('%s',fname_in));
disp(sprintf(' '));
disp(sprintf(' Now that the code has finished, we can load the output-ordering: '));
tmp_xdrop_fname = sprintf('%s%s_%s_out_xdrop.txt',path_base,prefix_base,frwd_vs_back);
disp(sprintf(' stored as %s',tmp_xdrop_fname));
out_xdrop = textread(tmp_xdrop_fname);
disp(sprintf(' This output-file lists the row and column indices in the order that they were eliminated. '));
disp(sprintf(' Each line of text within the output-file has two numbers, the first corresponding to  '));
disp(sprintf(' the row-index eliminated at that iteration of the algorithm (with a ''-1'' implying no elimination), '));
disp(sprintf(' and the second corresponding to the column-index eliminated at that iteration of the '));
disp(sprintf(' algorithm (again, with ''-1'' implying that no column was eliminated at that iteration). '));
disp(sprintf(' The convention we have chosen here is to list the eliminated row-indices using numerical '));
disp(sprintf(' values from 0 to (# of patients)-1. '));
disp(sprintf(' Similarly, we choose to represent the eliminated col-indices from 0 to (# of genes)-1. '));
disp(sprintf(' Note that some rows and columns are never eliminated, and are thus not listed in the '));
disp(sprintf(' output-file. This can easily be corrected for by reading all the eliminated indices, and '));
disp(sprintf(' then adding back in the unlisted indices; which we do now. '));
rdrop = out_xdrop(:,1); cdrop = out_xdrop(:,2); 
rij = rdrop(find(rdrop>-1)); cij = cdrop(find(cdrop>-1));
rij = [rij(:) ; setdiff(find(Ainc_n_rind)-1,rij(:))]; rij = rij(end:-1:1);
cij = [cij(:) ; setdiff(g_ij(:)-1,cij(:))]; cij = cij(end:-1:1);

disp(sprintf(' '));
disp(sprintf(' After reading in the output-file, we can reorganize the original matrix based on this output-file. '));
disp(sprintf(' Specifically, we can order the rows based on the reverse-order in which they were eliminated,'));
disp(sprintf(' and do the same thing for the columns. '));
disp(sprintf(' This allows us to organize the original case-matrix so that the longest-retained row- and '));
disp(sprintf(' column-indices are placed in the upper left corner. '));
disp(sprintf(' The control-matrix is shown below the case-matrix, with column-indices rearranged. '));

figure;
pgap = 32;ggap = 128;
cov_tmp = zeros(size(cov_mat,1),size(cov_mat,2)*ggap); for nc=1:size(cov_mat,2); cov_tmp(:,(nc-1)*ggap+(1:ggap)) = repmat(cov_mat(:,nc),1,ggap); end;%for nc=1:size(cov_mat,2);
cmap_autumn = colormap('autumn'); cmap_autumn(1,:) = [1,1,1]; cmap_autumn(end-1,:) = 0.85*[1,1,1]; cmap_autumn(end,:) = 0.00*[1,1,1]; clen = size(cmap_autumn,1); colormap(cmap_autumn);
cbot=-1.5; ctop=1.5; cvals = linspace(cbot,ctop,clen); cbotl=cvals(2); ctopl = cvals(end-2); cdiff = cvals(end)-cvals(end-1); 
Dtmp = Aall_n(1+rij,1+cij) ; cov_tmp_D = cov_tmp(1+rij,:);
Xtmp = Aall_n( x_ij,1+cij) ; cov_tmp_X = cov_tmp( x_ij,:);
subplot(1,2,1);
imagesc([...
  min(ctopl,max(cbotl,Dtmp)) , cbot*ones(length(1+rij),ggap) , ((ctopl+cdiff+cdiff*cov_tmp_D)) ; ...
  cbot*ones(pgap,ngenes) , cbot*ones(pgap,ggap) , cbot*ones(pgap,size(cov_tmp,2)); ...
  min(ctopl,max(cbotl,Xtmp)) , cbot*ones(length( x_ij),ggap) , ((ctopl+cdiff+cdiff*cov_tmp_X)) ; ...
	 ],[cbot,ctop]);
set(gca,'Xtick',[],'Ytick',[]); title('retained rows/cols (upper left) organized via output of lakcluster');
pgap = 32;ggap = 16;
cov_tmp = zeros(size(cov_mat,1),size(cov_mat,2)*ggap); for nc=1:size(cov_mat,2); cov_tmp(:,(nc-1)*ggap+(1:ggap)) = repmat(cov_mat(:,nc),1,ggap); end;%for nc=1:size(cov_mat,2);
Dtmp = Aall_n(1+rij(1:64),1+cij(1:196)) ; cov_tmp_D = cov_tmp(1+rij(1:64),:);
subplot(1,2,2);
imagesc([...
  min(ctopl,max(cbotl,Dtmp)) , cbot*ones(length(1+rij(1:64)),ggap) , ((ctopl+cdiff+cdiff*cov_tmp_D)) ; ...
	 ],[cbot,ctop]);
set(gca,'Xtick',[],'Ytick',[]); title('zoom into upper 64-x-196 corner');

disp(sprintf(' '));
disp(sprintf(' Paused (push a key). '));
pause();
disp(sprintf(' '));
disp(sprintf(' We can clearly see some sort of low-rank bicluster in the upper left corner, implying that the '));
disp(sprintf(' tutorial_lakcluster_1 algorithm found some kind of bicluster. But which one?'));
disp(sprintf(' In order to show which bicluster was retained, we plot, as a function of iteration, the '));
disp(sprintf(' overlap fraction of the retained rows (solid) and columns (dashed) with the various biclusters: '));
disp(sprintf(' bicluster 0: red '));
disp(sprintf(' bicluster 1: green '));
disp(sprintf(' bicluster 2: green '));
disp(sprintf(' bicluster 3: black '));
disp(sprintf(' '));
disp(sprintf(' Note that, as the algorithm progresses, the rows/columns which are retained (left side) strongly '));
disp(sprintf(' overlap with bicluster 0 (red), whereas bicluster 3 (black) is mostly eliminated fairly early on. '));

figure;
nbins = 128;
f_ra = linspace(0,1,nbins);
r_ra = max(1,min(length(d_ij),1+floor(f_ra.*length(d_ij))));
c_ra = max(1,min(ngenes,1+floor(f_ra.*ngenes)));
f_C0_r = zeros(nbins,1); f_C0_c = zeros(nbins,1);
f_C1_r = zeros(nbins,1); f_C1_c = zeros(nbins,1);
f_C2_r = zeros(nbins,1); f_C2_c = zeros(nbins,1);
f_C3_r = zeros(nbins,1); f_C3_c = zeros(nbins,1);
for nl=1:length(f_ra);
tmp_r = 1+rij(1:r_ra(nl)); tmp_c = 1+cij(1:c_ra(nl));
f_C0_r(nl) = length(intersect(tmp_r,d_ij(p_d_ij(1:C0rows))))/C0rows;
f_C0_c(nl) = length(intersect(tmp_c,g_ij(p0_g_ij(1:C0cols))))/C0cols;
f_C1_r(nl) = length(intersect(tmp_r,d_ij(c1d_ord(1:C1rows))))/C1rows;
f_C1_c(nl) = length(intersect(tmp_c,g_ij(p1_g_ij(1:C1cols))))/C1cols;
f_C2_r(nl) = length(intersect(tmp_r,d_ij(c2d_ord(1:C2rows))))/C2rows;
f_C2_c(nl) = length(intersect(tmp_c,g_ij(p2_g_ij(1:C2cols))))/C2cols;
f_C3_r(nl) = length(intersect(tmp_r,d_ij(q_d_ij(1:C3rows/2))))/(C3rows/2);
f_C3_c(nl) = length(intersect(tmp_c,g_ij(p3_g_ij(1:C3cols))))/C3cols;
end;%for nl=1:length(f_ra);
subplot(1,1,1);
hold on;
plot(f_ra,f_C0_r,'r-',f_ra,f_C1_r,'g-',f_ra,f_C2_r,'g-',f_ra,f_C3_r,'k-');
plot(f_ra,f_C0_c,'r--',f_ra,f_C1_c,'g--',f_ra,f_C2_c,'g--',f_ra,f_C3_c,'k--');
hold off;
xlabel('retained <--> discarded');
ylabel('fraction overlap');
legend('bc0','bc1','bc2','bc3','Location','SouthEast');
title('fraction of each bicluster retained (row=solid, column=dashed)');

disp(sprintf(' '));
disp(sprintf(' Paused (push a key). '));
pause();
disp(sprintf(' '));
disp(sprintf(' Now we quantify the output of the algorithm in a manner that does not rely on prior knowledge '));
disp(sprintf(' of the embedded biclusters. In order to do this we examine the iterations taken by the algorithm, '));
disp(sprintf(' each of which involves the calculation of ''row-scores'' and ''column-scores''. These scores '));
disp(sprintf(' estimate how likely it is that each row/column participates in a bicluster. '));
disp(sprintf(' These scores are corrected for covariates, to ensure that they emphasize biclusters which are '));
disp(sprintf(' balanced with respect to the covariates. '));
disp(sprintf(' As the algorithm proceeds, rows and columns with low scores are eliminated, and those with high '));
disp(sprintf(' scores are retained. '));
disp(sprintf(' Thus, to quantify how well the algorithm is proceeding, we could, at each iteration, measure the '));
disp(sprintf(' average scores of the remaining rows and columns. '));
disp(sprintf(' If these average scores are high, then the algorithm is focusing on a (covariate-balanced) bicluster; '));
disp(sprintf(' on the other hand, if these average scores are low, then the algorithm has not yet found anything. '));
disp(sprintf(' '));
tmp_trace_fname = sprintf('%s%s_%s_out_trace.txt',path_base,prefix_base,frwd_vs_back);
disp(sprintf(' These average row- and column-scores have actually already been calculated, and are stored within the '));
disp(sprintf(' output-file %s',tmp_trace_fname));
disp(sprintf(' This output-trace contains 1 line per iteration of the algorithm, with the numbers on that line listing:'));
disp(sprintf(' first number: the iteration number. '));
disp(sprintf(' second number: the number of case-patients remaining. '));
disp(sprintf(' third number: the number of genes remaining. '));
disp(sprintf(' fourth number: the average row-score across the remaining patients+genes. '));
disp(sprintf(' fifth number: the average column-score across the remaining patient+genes. '));
disp(sprintf(' sixth number: the number of covariate-categories remaining across the patients. '));
disp(sprintf(' '));
disp(sprintf(' As mentioned above, this output-trace can be used to assess whether or not the algorithm has succeeded. '));
disp(sprintf(' Indeed, in the case of real gene-expression data, we will use these output-trace values to estimate '));
disp(sprintf(' a p-value for the output-ordering (i.e., a p-value for any bicluster found). '));
disp(sprintf(' '));
disp(sprintf(' We plot this trace in the following figure:'));
[tmp] = textread(tmp_trace_fname);
figure;plot(tmp(:,1),tmp(:,4),'r.-',tmp(:,1),tmp(:,5),'b.-'); xlim([tmp(1,1),tmp(end,1)]); ylim([0,1.1]); legend('row','col','Location','SouthEast'); title('average score'); ylabel('average score'); xlabel('iteration');
disp(sprintf(' '));
disp(sprintf(' Paused (push a key). '));
pause();
disp(sprintf(' '));
disp(sprintf(' However, for the purposes of this tutorial we go a step further: '));
disp(sprintf(' we first use the output of the algorithm to reorder the data (as shown in an earlier '));
disp(sprintf(' figure), and then measure the average scores of each m-by-n ''corner submatrix'' C(m,n), formed '));
disp(sprintf(' from the m longest-retained rows and the n longest-retained columns. '));
disp(sprintf(' '));
disp(sprintf(' We embark on this task here by calling tutorial_lakcluster_1 with the loopS_flag on. '));
disp(sprintf(' For expedience we calculate these average scores for only the small corner submatrices, '));
disp(sprintf(' (of size < 128-by-512) rather than for all of them.'));
disp(sprintf(' Specifically, given the output-ordering of tutorial_lakcluster_1, we consider, for each (m,n), '));
disp(sprintf(' the corner-submatrix C(m,n) (comprising the m longest-retained rows and n longest-retained columns). '));
disp(sprintf(' Given C(m,n), we calculate: '));
disp(sprintf(' LR2_ra(m,n): which is the average row-score calculated across the m rows of C(m,n), and '));
disp(sprintf(' LC2_ra(m,n): which is the average column-score calculated across the n column of C(m,n), and '));
disp(sprintf(' Lx2_ra(m,n): which is the number of covariate-categories remaining within the m rows of C(m,n). '));
disp(sprintf(' '));
tmp_loopS_fname = sprintf('%s%s_%s_out_loopS.mat',path_base,prefix_base,frwd_vs_back);
disp(sprintf(' This file is stored as: %s',tmp_loopS_fname));
load_flag=0; if exist(tmp_loopS_fname,'file'); load_flag = input(' File exists: reload? [default = 1]: '); if isempty(load_flag); load_flag==1; end; end;
if ~load_flag; tutorial_lakcluster_1(sprintf('%s',fname_in),1,128,196); end;
load(tmp_loopS_fname);
disp(sprintf(' '));
disp(sprintf(' These arrays are shown in the following figure, with the first subplot showing LR2, the second shown LC2 '));
disp(sprintf(' (with the number of remaining covariate-categories plotted as a function of m to the far right. '));
disp(sprintf(' Note that the average scores are quite high for corner submatrices that are around 60-by-150. '));
disp(sprintf(' This protruding shape in the LR2 and LC2 arrays is indicative of revealed structure. '));

figure; cla;
subplot(1,5,[1,2]); imagesc(LR2_ra,[0,1]); colorbar; title('average row score'); set(gca,'Ytick',[1,30,60,90],'Yticklabel',[1,30,60,90],'Xtick',[1,50,100],'Xticklabel',[1,50,100]); xlabel('columns n'); ylabel('rows m');
subplot(1,5,[3,4]); imagesc(LC2_ra,[0,1]); colorbar; title('average col score'); set(gca,'Ytick',[1,30,60,90],'Yticklabel',[1,30,60,90],'Xtick',[1,50,100],'Xticklabel',[1,50,100]); xlabel('columns n'); ylabel('rows m');
subplot(1,5,5); stairs(Lx2_ra(end:-1:1,1),1:127);xlim([0,5]);ylim([0,127]);set(gca,'Ytick',128-[90,60,30,1],'Yticklabel',[90,60,30,1],'Xtick',[1,2,3,4],'Xticklabel',[1,2,3,4]); title('# covariates remaining');
xlabel('covariate-categories remaining'); ylabel('rows m');

disp(sprintf(' '));
disp(sprintf(' Paused (push a key). '));
pause();
disp(sprintf(' '));
disp(sprintf(' If we wanted to, we could threshold these average-scores and extract something like this 60-by-150 '));
disp(sprintf(' submatrix, labeling it as a bicluster. '));
disp(sprintf(' To do this systematically, we run the program tutorial_loopS_threshold.m. '));
disp(sprintf(' This program requires as input an array of average-scores to threshold. '));
disp(sprintf(' We choose to use LR2_ra because, as in the case of gene-expression analysis, we often have fewer rows '));
disp(sprintf(' (i.e., patients) than columns (i.e., genes); the row-scores are typically more informative than the '));
disp(sprintf(' column-scores. '));
disp(sprintf(' Given these average-scores, we are looking for a corner matrix that scores highly, and is ''big''. '));
disp(sprintf(' To quantify ''bigness'', we search for the corner-submatrix that has the largest possible area, '));
disp(sprintf(' subject to the constraint that the average score of that corner-submatrix is at least some minimum '));
disp(sprintf(' value (denoted by LR2_thr). '));
disp(sprintf(' We define LR2_thr to be 70%% of the way from the minimum to the maximum value of LR2_ra. '));
disp(sprintf(' In addition to this ''magic number'' of 70%% (denoted rthr -- for relative threshold), '));
disp(sprintf(' the program tutorial_loopS_threshold.m also requires a few other parameters: '));
disp(sprintf(' Lx2_ra: only corner-submatrices with the full number of covariate-categories will be considered '));
disp(sprintf(' ar_min, ar_max: only corner-submatrices with relative-aspect-ratios between ar_min and ar_max will be considered '));
disp(sprintf(' plot_flag: if this is set to 1 then we plot the output, as we do here. '));
disp(sprintf(' '));
rthr = 0.7; ar_min = 0.25 ; ar_max = 1.0/ar_min ;
[output,ik_out,ij_out] = tutorial_loopS_threshold(LR2_ra,Lx2_ra,rthr,ar_min,ar_max,1);
disp(sprintf(' '));
disp(sprintf(' Within this plot the valid corners (i.e., within the aspect-ratio bounds and above the threshold) '));
disp(sprintf(' are shown in white, with the maximal area corner marked by an ''x''. '));
disp(sprintf(' This maximal-area corner is 54-by-159. '));
disp(sprintf(' Note that the colorbar for this plot (LR2_min to LR2_max) is different than the previous figure (0,1). '));
disp(sprintf(' '));
disp(sprintf(' Paused (push a key). '));
pause();
disp(sprintf(' '));
disp(sprintf(' At this point we can extract this 54-by-159 submatrix and plot it in relation to the controls: '));
disp(sprintf(' As one can see, this bicluster consists of genes that are strongly co-expressed '));
disp(sprintf(' across the extracted subset of patients, relative to the full set of controls '));
disp(sprintf(' (to make this point clear we have rearranged the rows and columns of each matrix). '));
disp(sprintf(' In the context of gene-expression analysis we may want to remove this submatrix from the data '));
disp(sprintf(' (e.g., by replacing it with random values drawn from the rest of the matrix) and then search '));
disp(sprintf(' for a second bicluster. This is carried out towards the end of tutorial_w1.m by passing '));
disp(sprintf(' additional arguments to tutorial_lakcluster_1.m via an expanded input file. '));

if (output>0);
ncij = 1:nc_ra(ik_out); nrij = 1:nr_ra(ij_out);
rij_sub = rij(nrij); cij_sub = cij(ncij);
pgap = 32;ggap = 4;
cov_tmp = zeros(size(cov_mat,1),size(cov_mat,2)*ggap); for nc=1:size(cov_mat,2); cov_tmp(:,(nc-1)*ggap+(1:ggap)) = repmat(cov_mat(:,nc),1,ggap); end;%for nc=1:size(cov_mat,2);
cmap_autumn = colormap('autumn'); cmap_autumn(1,:) = [1,1,1]; cmap_autumn(end-1,:) = 0.85*[1,1,1]; cmap_autumn(end,:) = 0.00*[1,1,1]; clen = size(cmap_autumn,1); colormap(cmap_autumn);
cbot=-1.5; ctop=1.5; cvals = linspace(cbot,ctop,clen); cbotl=cvals(2); ctopl = cvals(end-2); cdiff = cvals(end)-cvals(end-1); 
Dtmp = Aall_n(1+rij_sub,1+cij_sub) ; cov_tmp_D = cov_tmp(1+rij_sub,:);
Xtmp = Aall_n( x_ij,1+cij_sub) ; cov_tmp_X = cov_tmp( x_ij,:);
[prA,pc,prX] = tutorial_hcluster(Dtmp,Xtmp); Dtmp = Dtmp(prA,pc); Xtmp = Xtmp(prX,pc); cov_tmp_D = cov_tmp_D(prA,:); cov_tmp_X = cov_tmp_X(prX,:);
imagesc([...
  min(ctopl,max(cbotl,Dtmp)) , cbot*ones(length(1+rij_sub),ggap) , ((ctopl+cdiff+cdiff*cov_tmp_D)) ; ...
  cbot*ones(pgap,length(cij_sub)) , cbot*ones(pgap,ggap) , cbot*ones(pgap,size(cov_tmp,2)); ...
  min(ctopl,max(cbotl,Xtmp)) , cbot*ones(length( x_ij),ggap) , ((ctopl+cdiff+cdiff*cov_tmp_X)) ; ...
	 ],[cbot,ctop]);
set(gca,'Xtick',[],'Ytick',[]); title(sprintf('extracted submatrix (top: corresponding to relative threshold %0.2f) vs controls (bottom)',rthr));
end;%if (output>0);
