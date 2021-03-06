function tutorial_summarize(path_base,gen_fname,frwd_vs_back);
% designed to summarize data generated by helper_collate_6.m;
% loads bc, xdrop and trace ;
% try with: ;
% tutorial_summarize('/data/rangan/dir_bcc/dir_tutorial_081915/dir_GSE17536/','GSE17536_n0x','frwd');

disp_flag=1;

bitj=16; 
prefix_base = sprintf('%s',gen_fname); prefix = prefix_base; 
path_pre_prm = sprintf('%sdir_%s_prm/',path_base,prefix_base);
path_FIGS = sprintf('%sdir_FIGS/',path_base); if ~exist(path_FIGS,'dir'); mkdir(path_base,'dir_FIGS'); end;
path_use = path_base;
prefix = prefix_base; 
tutorial_w1_setup;

load(sprintf('%s%s_%s_collate.mat',path_FIGS,prefix_base,frwd_vs_back),'Omax','ox_ra','O_ra','NO','Pmax','ex_ra_','path_pre_prm','P_xxxx_ra_','NP_xxxx_','bc_size','tr_base_ra','tr_ra_','tc_base_ra','tc_ra_','nc_base_ra','nc_ra_','cij_base_ra','tr_avg_pval_','tr_zmx_pval_','tr_max_pval_','tr_avg_','tr_zmx_','tr_max_','tr_base_avg_ra','tr_base_zmx_ra','tr_base_max_ra','tr_cmb_pval_');

disp_flag=1;
fname_out = sprintf('%sdir_txt/%s_%s_summarize.txt',path_use,prefix,frwd_vs_back);
fid_out = fopen(fname_out,'w');
fname_sum = sprintf('%s../BCsum.m',path_use);
gse_tab = ...
  ~isempty(strfind(path_use,'GSE17536'))* 2+ ...
  0;
n_tab = ...
~isempty(strfind(gen_fname,'_n0__'))* 0 + ...
~isempty(strfind(gen_fname,'_n0x_'))* 0 + ...
~isempty(strfind(gen_fname,'_n1__'))* 1 + ...
~isempty(strfind(gen_fname,'_n1x_'))* 1 + ...
~isempty(strfind(gen_fname,'_n2__'))* 2 + ...
~isempty(strfind(gen_fname,'_n2x_'))* 2 + ...
~isempty(strfind(gen_fname,'_nA__'))* 4 + ...
~isempty(strfind(gen_fname,'_nAx_'))* 4 + ...
  0;
x_tab = ...
~isempty(strfind(gen_fname,'_n0__'))* 0 + ...
~isempty(strfind(gen_fname,'_n0x_'))* 1 + ...
~isempty(strfind(gen_fname,'_n1__'))* 0 + ...
~isempty(strfind(gen_fname,'_n1x_'))* 1 + ...
~isempty(strfind(gen_fname,'_n2__'))* 0 + ...
~isempty(strfind(gen_fname,'_n2x_'))* 1 + ...
~isempty(strfind(gen_fname,'_nA__'))* 0 + ...
~isempty(strfind(gen_fname,'_nAx_'))* 1 + ...
  0;
fb_tab = ~isempty(strfind(frwd_vs_back,'frwd'));
fid_sum = fopen(fname_sum,'a');
for no=1:NO;
p_tmp = tr_cmb_pval_(no);
tmpstr = sprintf('%% bc%d [%d-x-%d]: p_val: %0.4f',no-1,bc_size(1,no),bc_size(2,no),p_tmp);
fprintf(fid_sum,'BCsum(%d,%d,%d,%d,%d,1:4) = [%d,%d,0,0];\n',1+gse_tab,1+n_tab,1+x_tab,1+fb_tab,1+O_ra(no)-1,bc_size(1,no),bc_size(2,no));
fprintf(fid_sum,'BCsum(%d,%d,%d,%d,%d,4 + (%d-1)*8 + (1:8)) = [%f,%f,%f,%f,%f,%f,%f,%f];\n',1+gse_tab,1+n_tab,1+x_tab,1+fb_tab,1+O_ra(no)-1,0,[tr_max_pval_(no),tr_zmx_pval_(no),tr_avg_pval_(no)],p_tmp,0);
if disp_flag;
if 1 | p_tmp(1)<0.1 ;
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(fid_out,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
disp(tmpstr);
fprintf(fid_out,'%s\n',tmpstr);
no_flag = O_ra(no)-1;
if no_flag==0; tmp_genri_fname = sprintf('%sdir_txt/%s_%s_genri.txt',path_use,prefix,frwd_vs_back);
 else; tmp_genri_fname = sprintf('%sdir_txt/%s_%s_%d_genri.txt',path_use,prefix,frwd_vs_back,no_flag); end;
disp('%genri---------------------------------------------------------%');
fprintf(fid_out,'%%genri---------------------------------------------------------%%\n');
fid_tmp = fopen(tmp_genri_fname,'r'); for tmp_ij=1:10;tmpln=fgetl(fid_tmp); if ischar(tmpln); disp(tmpln); fprintf(fid_out,sprintf('%s\n',tmpln)); end; end; %for tmp_ij=1:10;fclose(fid_tmp);
fid_tmp = fopen(tmp_genri_fname,'r'); tmp_C = textscan(fid_tmp,'%s %f %f %d %d','headerlines',1); fclose(fid_tmp); 
if (length(tmp_C{2})>0); 
if (tmp_C{2}(1)<0.001); fprintf(fid_sum,'BCsum(%d,%d,%d,%d,%d,3) = [2];\n',1+gse_tab,1+n_tab,1+x_tab,1+fb_tab,1+O_ra(no)-1);
 else fprintf(fid_sum,'BCsum(%d,%d,%d,%d,%d,3) = [1];\n',1+gse_tab,1+n_tab,1+x_tab,1+fb_tab,1+O_ra(no)-1); end;
end;%if (length(tmp_C{2})>0); 
if no_flag==0; tmp_gslim_fname = sprintf('%sdir_txt/%s_%s_gslim.txt',path_use,prefix,frwd_vs_back);
 else; tmp_gslim_fname = sprintf('%sdir_txt/%s_%s_%d_gslim.txt',path_use,prefix,frwd_vs_back,no_flag); end;
disp('%gslim---------------------------------------------------------%');
fprintf(fid_out,'%%gslim---------------------------------------------------------%%\n');
fid_tmp = fopen(tmp_gslim_fname,'r'); for tmp_ij=1:10;tmpln=fgetl(fid_tmp); if ischar(tmpln); disp(tmpln); fprintf(fid_out,sprintf('%s\n',tmpln)); end; end; %for tmp_ij=1:10;fclose(fid_tmp);
fid_tmp = fopen(tmp_gslim_fname,'r'); tmp_C = textscan(fid_tmp,'%s %f %f %d %d','headerlines',1); fclose(fid_tmp); 
if (length(tmp_C{2})>0); 
if (tmp_C{2}(1)<0.001); fprintf(fid_sum,'BCsum(%d,%d,%d,%d,%d,4) = [2];\n',1+gse_tab,1+n_tab,1+x_tab,1+fb_tab,1+O_ra(no)-1);
 else fprintf(fid_sum,'BCsum(%d,%d,%d,%d,%d,4) = [1];\n',1+gse_tab,1+n_tab,1+x_tab,1+fb_tab,1+O_ra(no)-1); end;
end;%if (length(tmp_C{2})>0); 
end;%if 1 | p_tmp(1)<0.1 ;
end;%if disp_flag;
end;%for no=1:NO;
if disp_flag;
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(fid_out,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
end;%if disp_flag;
fclose(fid_sum);
fclose(fid_out);
