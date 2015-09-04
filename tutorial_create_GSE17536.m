clear all;
fpre = 'GSE17536'; ppre = 'GPL570'; path_put = sprintf('./dir_%s/',fpre); fpre_put = sprintf('%s%s',path_put,fpre);
if (~exist(path_put,'dir')); command_string = sprintf('mkdir %s;',path_put); system(command_string); end;
if (exist(sprintf('./%s.%s.pcl',fpre,ppre),'file') & ~exist(sprintf('%s.%s.pcl',fpre_put,ppre),'file')); command_string = sprintf('mv ./%s.%s.pcl %s.%s.pcl ;',fpre,ppre,fpre_put,ppre); system(command_string); end;
if (exist(sprintf('./%s.gsms.txt',fpre),'file') & ~exist(sprintf('%s.gsms.txt',fpre_put),'file')); command_string = sprintf('mv ./%s.gsms.txt %s.gsms.txt ;',fpre,fpre_put); system(command_string); end;

disp(sprintf(' This file reads from the files: '));
disp(sprintf(' %s.gsms.txt and %s.%s.pcl, ',fpre,fpre,ppre));
disp(sprintf(' which hold the annotations and data for the %s series, respectively. ',fpre));

disp(sprintf('  '));
disp(sprintf(' First we use the linux shell to extract the following information from the data files: '));
disp(sprintf(' %s_dat_order.txt: header of %s.%s.pcl; contains patient indices within pcl file ',fpre_put,fpre_put,ppre));
disp(sprintf(' %s_hdr_order.txt: header of %s.gsms.txt; contains patient indices within gsmsm file ',fpre_put,fpre_put));
disp(sprintf(' %s___male_index.txt; male patients ',fpre_put));
disp(sprintf(' %s_female_index.txt; female patients ',fpre_put));
disp(sprintf(' %s_ethn_b_index.txt; ethnicity black ',fpre_put));
disp(sprintf(' %s_ethn_c_index.txt; ethnicity caucasian ',fpre_put));
disp(sprintf(' %s_ethn_h_index.txt; ethnicity hispanic ',fpre_put));
disp(sprintf(' %s_ethn_o_index.txt; ethnicity other ',fpre_put));
disp(sprintf(' %s_over_d_index.txt; overall event --> death ',fpre_put));
disp(sprintf(' %s_over_n_index.txt; overall event --> no death ',fpre_put));
disp(sprintf(' %s_dss_de_index.txt; (disease specific survival; death from cancer): death ',fpre_put));
disp(sprintf(' %s_dss_no_index.txt; (disease specific survival; death from cancer): no death ',fpre_put));
disp(sprintf(' %s_dfs_NA_index.txt; (disease free survival; cancer recurrence): NA ',fpre_put));
disp(sprintf(' %s_dfs_no_index.txt; (disease free survival; cancer recurrence): no recurrence ',fpre_put));
disp(sprintf(' %s_dfs_re_index.txt; (disease free survival; cancer recurrence): recurrence ',fpre_put));
disp(sprintf(' %s____age_vals.txt; patient age ',fpre_put));
disp(sprintf(' %s___ajcc_vals.txt; ajcc_stage ',fpre_put));
disp(sprintf(' %s__grade_vals.txt; grade value ',fpre_put));

disp(sprintf('  '));
disp(sprintf(' We will later use the disease specific survival (death vs no death) to categorize the patients into ''cases'' and ''controls''. '));
disp(sprintf(' We will eventually search for biclusters within the cases (i.e., biclusters within patients that died of cancer). '));

linux_flag=1;%linux_flag=input('extract indices from headers on linux? (1 for yes, 0 for no): ');
if linux_flag;
command = sprintf('head -1 %s.%s.pcl >! %s_dat_order.txt',fpre_put,ppre,fpre_put); system(command);
command = sprintf('cut -f1 -d''|'' %s.gsms.txt >! %s_hdr_order.txt',fpre_put,fpre_put); system(command);
command = sprintf('grep -n                                                       ''|gender: female|'' %s.gsms.txt | cut -f1 -d: >! %s_female_index.txt',fpre_put,fpre_put); system(command); 
command = sprintf('grep -n                                                         ''|gender: male|'' %s.gsms.txt | cut -f1 -d: >! %s___male_index.txt',fpre_put,fpre_put); system(command); 
command = sprintf('grep -n                                                       ''ethnicity: other'' %s.gsms.txt | cut -f1 -d: >! %s_ethn_o_index.txt',fpre_put,fpre_put); system(command); 
command = sprintf('grep -n                                                       ''ethnicity: black'' %s.gsms.txt | cut -f1 -d: >! %s_ethn_b_index.txt',fpre_put,fpre_put); system(command); 
command = sprintf('grep -n                                                    ''ethnicity: hispanic'' %s.gsms.txt | cut -f1 -d: >! %s_ethn_h_index.txt',fpre_put,fpre_put); system(command); 
command = sprintf('grep -n                                                   ''ethnicity: caucasian'' %s.gsms.txt | cut -f1 -d: >! %s_ethn_c_index.txt',fpre_put,fpre_put); system(command); 
command = sprintf('grep -n                          ''|overall_event (death from any cause): death|'' %s.gsms.txt | cut -f1 -d: >! %s_over_d_index.txt',fpre_put,fpre_put); system(command); 
command = sprintf('grep -n                       ''|overall_event (death from any cause): no death|'' %s.gsms.txt | cut -f1 -d: >! %s_over_n_index.txt',fpre_put,fpre_put); system(command); 
command = sprintf('grep -n      ''|dss_event (disease specific survival; death from cancer): death|'' %s.gsms.txt | cut -f1 -d: >! %s_dss_de_index.txt',fpre_put,fpre_put); system(command); 
command = sprintf('grep -n   ''|dss_event (disease specific survival; death from cancer): no death|'' %s.gsms.txt | cut -f1 -d: >! %s_dss_no_index.txt',fpre_put,fpre_put); system(command); 
command = sprintf('grep -n  ''|dfs_event (disease free survival; cancer recurrence): no recurrence|'' %s.gsms.txt | cut -f1 -d: >! %s_dfs_no_index.txt',fpre_put,fpre_put); system(command); 
command = sprintf('grep -n     ''|dfs_event (disease free survival; cancer recurrence): recurrence|'' %s.gsms.txt | cut -f1 -d: >! %s_dfs_re_index.txt',fpre_put,fpre_put); system(command); 
command = sprintf('grep -n             ''|dfs_event (disease free survival; cancer recurrence): NA|'' %s.gsms.txt | cut -f1 -d: >! %s_dfs_NA_index.txt',fpre_put,fpre_put); system(command); 
command = sprintf('cat %s.gsms.txt | cut -f4 -d''|'' | cut -f2 -d'':''            >! %s____age_vals.txt',fpre_put,fpre_put); system(command);
command = sprintf('cat %s.gsms.txt | cut -f7 -d''|'' | cut -f2 -d'':''            >! %s___ajcc_vals.txt',fpre_put,fpre_put); system(command);
command = sprintf('cat %s.gsms.txt | cut -f8 -d''|'' | cut -f2 -d'':'' | cut -c-2 >! %s__grade_vals.txt',fpre_put,fpre_put); system(command);
%return;
end;%if linux_flag;

disp(sprintf('  '));
disp(sprintf(' Now we read the patient-ids from the dat_order and hdr_order files, noting that the order may not be the same!  '));
disp(sprintf(' Moreover, there are more patient-ids (177) in the dat_order file than in the hdr_order file (175), '));
disp(sprintf(' implying that some patients were not provided with annotations. '));
disp(sprintf(' In this case the ''unannotated'' patients are GSM437093 and GSM437094; the first 2 in the dat_order list. '));
fid = fopen(sprintf('%s_dat_order.txt',fpre_put)); dat_order = textscan(fid,'%s',177); fclose(fid); dat_order = dat_order{1}; for nl=1:length(dat_order); tmp=dat_order(nl); dat_val(nl)=str2num(tmp{1}(4:end)); end;
fid = fopen(sprintf('%s_hdr_order.txt',fpre_put)); hdr_order = textscan(fid,'%s',175); fclose(fid); hdr_order = hdr_order{1}; for nl=1:length(hdr_order); tmp=hdr_order(nl); hdr_val(nl)=str2num(tmp{1}(4:end)); end;
[iv,hdr_ij,dat_ij] = intersect(hdr_val,dat_val,'stable');
npats = length(hdr_ij);
cov___male_ij = textread(sprintf('%s___male_index.txt',fpre_put));cov___male = zeros(npats,1); cov___male(cov___male_ij) = 1;
cov_female_ij = textread(sprintf('%s_female_index.txt',fpre_put));cov_female = zeros(npats,1); cov_female(cov_female_ij) = 1;
cov_ethn_o_ij = textread(sprintf('%s_ethn_o_index.txt',fpre_put));cov_ethn_o = zeros(npats,1); cov_ethn_o(cov_ethn_o_ij) = 1;
cov_ethn_b_ij = textread(sprintf('%s_ethn_b_index.txt',fpre_put));cov_ethn_b = zeros(npats,1); cov_ethn_b(cov_ethn_b_ij) = 1;
cov_ethn_h_ij = textread(sprintf('%s_ethn_h_index.txt',fpre_put));cov_ethn_h = zeros(npats,1); cov_ethn_h(cov_ethn_h_ij) = 1;
cov_ethn_c_ij = textread(sprintf('%s_ethn_c_index.txt',fpre_put));cov_ethn_c = zeros(npats,1); cov_ethn_c(cov_ethn_c_ij) = 1;
cov_over_d_ij = textread(sprintf('%s_over_d_index.txt',fpre_put));cov_over_d = zeros(npats,1); cov_over_d(cov_over_d_ij) = 1;
cov_over_n_ij = textread(sprintf('%s_over_n_index.txt',fpre_put));cov_over_n = zeros(npats,1); cov_over_n(cov_over_n_ij) = 1;
cov_dss_de_ij = textread(sprintf('%s_dss_de_index.txt',fpre_put));cov_dss_de = zeros(npats,1); cov_dss_de(cov_dss_de_ij) = 1;
cov_dss_no_ij = textread(sprintf('%s_dss_no_index.txt',fpre_put));cov_dss_no = zeros(npats,1); cov_dss_no(cov_dss_no_ij) = 1;
cov_dfs_no_ij = textread(sprintf('%s_dfs_no_index.txt',fpre_put));cov_dfs_no = zeros(npats,1); cov_dfs_no(cov_dfs_no_ij) = 1;
cov_dfs_re_ij = textread(sprintf('%s_dfs_re_index.txt',fpre_put));cov_dfs_re = zeros(npats,1); cov_dfs_re(cov_dfs_re_ij) = 1;
cov_dfs_NA_ij = textread(sprintf('%s_dfs_NA_index.txt',fpre_put));cov_dfs_NA = zeros(npats,1); cov_dfs_NA(cov_dfs_NA_ij) = 1;
cov___age_vals = textread(sprintf('%s____age_vals.txt',fpre_put));
cov__ajcc_vals = textread(sprintf('%s___ajcc_vals.txt',fpre_put));
cov_grade_vals = textread(sprintf('%s__grade_vals.txt',fpre_put)); 
hdr_val = hdr_val(hdr_ij);
disp(sprintf('  '));
disp(sprintf(' Now we determine our final list of patient-ids (Pnames) using the ordering from the hdr_order file. '));
disp(sprintf(' We then create covariate vectors aligned with this list of patients. '));
Pnames = hdr_val;
cov___male = cov___male(hdr_ij);
cov_female = cov_female(hdr_ij);
cov_ethn_o = cov_ethn_o(hdr_ij);
cov_ethn_b = cov_ethn_b(hdr_ij);
cov_ethn_h = cov_ethn_h(hdr_ij);
cov_ethn_c = cov_ethn_c(hdr_ij);
cov_over_d = cov_over_d(hdr_ij);
cov_over_n = cov_over_n(hdr_ij);
cov_dss_de = cov_dss_de(hdr_ij);
cov_dss_no = cov_dss_no(hdr_ij);
cov_dfs_no = cov_dfs_no(hdr_ij);
cov_dfs_re = cov_dfs_re(hdr_ij);
cov_dfs_NA = cov_dfs_NA(hdr_ij);
cov___age_vals = cov___age_vals(hdr_ij);
cov__ajcc_vals = cov__ajcc_vals(hdr_ij);
cov_grade_vals = cov_grade_vals(hdr_ij);
cov_grade1 = cov_grade_vals==1;
cov_grade2 = cov_grade_vals==2;
cov_grade3 = cov_grade_vals==3;
cov_age_up = cov___age_vals>65; 

disp(sprintf('  '));
disp(sprintf(' Now we read from the pcl file, noting that, after skipping the first line (which we read earlier as dat_order), '));
disp(sprintf(' we need to read 177+1=178 values per line. This is because we need to read the gene-id (a number) followed by '));
disp(sprintf(' the 177 values associated with each of the 177 patient-ids stored in the pcl file. '));
disp(sprintf(' Note that after we read all these values, we will restrict our attention to those patient-ids within the intersection '));
disp(sprintf(' of the hdr_order and dat_order files; in this case we ignore the GSM437093 and GSM437094 patients (i.e., the first 2 columns). '));
disp(sprintf(' We will also later transpose our data (i.e., B) so that the rows correspond to patients and the columns to genes. '));
fid = fopen(sprintf('%s.%s.pcl',fpre_put,ppre));
A = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','headerlines',1);
fclose(fid);
ngenes = length(A{1});
B = zeros(ngenes,length(A));
for np=1:length(A);
B(:,np) = A{np};
end;%for np=1:length(A);
Gnames = B(:,1);
B = B(:,2:end);
B = B(:,dat_ij);
disp(sprintf('  '));
disp(sprintf(' We can see that the raw data (B) is unevenly distributed, with many large values. '));
disp(sprintf(' This can be seen, e.g., by comparing histograms of B and a normalized version called C. '));
disp(sprintf(' This C was created by taking the log of the data and then normalizing each gene. '));
C = transpose(log(B));
for ng = 1:ngenes;
C(:,ng) = (C(:,ng) - mean(C(:,ng)))/std(C(:,ng));
end;%for ng = 1:ngenes;
plot_flag = input(' Plot histograms? [default = 0]: '); if isempty(plot_flag);plot_flag=0;end; if plot_flag; subplot(2,1,1);hist(B(:),linspace(4,16,128));xlim([4,16]);title('B distribution'); subplot(2,1,2);hist(C(:),linspace(-5,5,128));xlim([-5.5,5.5]); title('C distribution'); end;
disp(sprintf(' We use this C for visualization, but not for our biclustering. '));

disp(sprintf('  '));
disp(sprintf(' At this point we need to decide which covariates to use. '));
disp(sprintf(' Two of the covariates -- namely overall-death and cancer-recurrence -- correlate very strongly '));
disp(sprintf(' with disease specific death (i.e., dss_de correlates with over_d and dfs_re). '));
disp(sprintf(' We will not control for these two covariates, as these are not expected to be confounders. '));
disp(sprintf(' However, we will make an attempt to control for some of the other covariates. '));
disp(sprintf(' This is somewhat tricky, because a first glance does not reveal any obvious structure '));
disp(sprintf(' linking many of the other covariates, such as age, to disease-specific-death.  '));
disp(sprintf(' This can be quantified via: [h,p] = ttest2(cov___age_vals(cov_dss_de_ij),cov___age_vals(cov_dss_no_ij)), '));
[h,p] = ttest2(cov___age_vals(cov_dss_de_ij),cov___age_vals(cov_dss_no_ij));
disp(sprintf(' which yields an insignificant P-value of %0.2f. ',p));
disp(sprintf(' At this point we can look a little more closely at the data... '));

plot_flag = input(' Plot data in 3-dimensional PCA-space; highlighting various covariates? [default = 0]: '); if isempty(plot_flag);plot_flag=0;end;
if (plot_flag);
disp(sprintf(' Note that these are 3-d plots, and can be rotated '));
% simple PCA-plot; note that these are 3-d plots, and can be rotated ;
% it does not seem as though gender, grade, dfs or ethn plays a huge role in patient clustering ;
[U,S,V] = svds(C,3);
figure; cla;
hold on;
plot3(U(find(cov___male),1),U(find(cov___male),2),U(find(cov___male),3),'r.','Markersize',25);
plot3(U(find(cov_female),1),U(find(cov_female),2),U(find(cov_female),3),'b.','Markersize',25);
hold off;
axis vis3d;
title('sort by gender (ma,fe)');
figure; cla;
hold on;
plot3(U(find(cov_grade1),1),U(find(cov_grade1),2),U(find(cov_grade1),3),'r.','Markersize',25);
plot3(U(find(cov_grade2),1),U(find(cov_grade2),2),U(find(cov_grade2),3),'g.','Markersize',25);
plot3(U(find(cov_grade3),1),U(find(cov_grade3),2),U(find(cov_grade3),3),'b.','Markersize',25);
hold off;
axis vis3d;
title('sort by grade (1,2,3)');
figure; cla;
hold on;
plot3(U(find(cov_dfs_no),1),U(find(cov_dfs_no),2),U(find(cov_dfs_no),3),'b.','Markersize',25);
plot3(U(find(cov_dfs_re),1),U(find(cov_dfs_re),2),U(find(cov_dfs_re),3),'r.','Markersize',25);
plot3(U(find(cov_dfs_NA),1),U(find(cov_dfs_NA),2),U(find(cov_dfs_NA),3),'k.','Markersize',25);
hold off;
axis vis3d;
title('sort by dfs (no,re,NA)');
figure; cla;
hold on;
plot3(U(find(cov_ethn_o),1),U(find(cov_ethn_o),2),U(find(cov_ethn_o),3),'k.','Markersize',25);
plot3(U(find(cov_ethn_b),1),U(find(cov_ethn_b),2),U(find(cov_ethn_b),3),'b.','Markersize',25);
plot3(U(find(cov_ethn_h),1),U(find(cov_ethn_h),2),U(find(cov_ethn_h),3),'g.','Markersize',25);
plot3(U(find(cov_ethn_c),1),U(find(cov_ethn_c),2),U(find(cov_ethn_c),3),'r.','Markersize',25);
hold off;
axis vis3d;
title('sort by ethn (o,b,h,c)');
end;%if (plot_flag);

plot_flag = input(' Plot heatmaps of gene-expression data sorted by the various covariates (see far right)? [default = 0]: '); if isempty(plot_flag);plot_flag=0;end;
if (plot_flag);
pgap = 4;ggap = 128;
cmap_autumn = colormap('autumn'); cmap_autumn(1,:) = [1,1,1]; clen = size(cmap_autumn,1); colormap(cmap_autumn);
cbot=-1.5; ctop=1.5; cvals = linspace(cbot,ctop,clen); cbotl=cvals(2); ctopl = cvals(end-1); 
cov_tmp = [repmat(cov___male,1,ggap) , repmat(cov_female,1,ggap) , repmat(cov_over_d,1,ggap) , repmat(cov_over_n,1,ggap) , repmat(cov_dss_de,1,ggap) , repmat(cov_dss_no,1,ggap) , repmat(cov_dfs_no,1,ggap) , repmat(cov_dfs_re,1,ggap) , repmat(cov_dfs_NA,1,ggap) , repmat(cov_grade1,1,ggap) , repmat(cov_grade2,1,ggap) , repmat(cov_grade3,1,ggap) , repmat(cov_ethn_o,1,ggap) , repmat(cov_ethn_b,1,ggap) , repmat(cov_ethn_h,1,ggap) , repmat(cov_ethn_c,1,ggap)];
colormap(cmap_autumn);
imagesc([...
	    min(ctopl,max(cbotl,C(find(cov___male),:))) , cbot*ones(length(find(cov___male)),ggap) , cov_tmp(find(cov___male),:) ; ...
	     cbot*ones(pgap,ngenes) , cbot*ones(pgap,ggap) , cbot*ones(pgap,size(cov_tmp,2)); ...
	    min(ctopl,max(cbotl,C(find(cov_female),:))) , cbot*ones(length(find(cov_female)),ggap) , cov_tmp(find(cov_female),:) ; ...
	 ],[cbot,ctop]); axis off;
title('sort by gender (ma,fe)');
figure;colormap(cmap_autumn);
imagesc([...
	    min(ctopl,max(cbotl,C(find(cov_grade1),:))) , cbot*ones(length(find(cov_grade1)),ggap) , cov_tmp(find(cov_grade1),:) ; ...
	     cbot*ones(pgap,ngenes) , cbot*ones(pgap,ggap) , cbot*ones(pgap,size(cov_tmp,2)); ...
	    min(ctopl,max(cbotl,C(find(cov_grade2),:))) , cbot*ones(length(find(cov_grade2)),ggap) , cov_tmp(find(cov_grade2),:) ; ...
	     cbot*ones(pgap,ngenes) , cbot*ones(pgap,ggap) , cbot*ones(pgap,size(cov_tmp,2)); ...
	    min(ctopl,max(cbotl,C(find(cov_grade3),:))) , cbot*ones(length(find(cov_grade3)),ggap) , cov_tmp(find(cov_grade3),:) ; ...
	 ],[cbot,ctop]); axis off;
title('sort by grade (1,2,3)');
figure;colormap(cmap_autumn);
imagesc([...
	    min(ctopl,max(cbotl,C(find(cov_dfs_no),:))) , cbot*ones(length(find(cov_dfs_no)),ggap) , cov_tmp(find(cov_dfs_no),:) ; ...
	     cbot*ones(pgap,ngenes) , cbot*ones(pgap,ggap) , cbot*ones(pgap,size(cov_tmp,2)); ...
	    min(ctopl,max(cbotl,C(find(cov_dfs_re),:))) , cbot*ones(length(find(cov_dfs_re)),ggap) , cov_tmp(find(cov_dfs_re),:) ; ...
	     cbot*ones(pgap,ngenes) , cbot*ones(pgap,ggap) , cbot*ones(pgap,size(cov_tmp,2)); ...
	    min(ctopl,max(cbotl,C(find(cov_dfs_NA),:))) , cbot*ones(length(find(cov_dfs_NA)),ggap) , cov_tmp(find(cov_dfs_NA),:) ; ...
	 ],[cbot,ctop]); axis off;
title('sort by dfs (no, re, NA)');
figure;colormap(cmap_autumn);
imagesc([...
	    min(ctopl,max(cbotl,C(find(cov_ethn_o),:))) , cbot*ones(length(find(cov_ethn_o)),ggap) , cov_tmp(find(cov_ethn_o),:) ; ...
	     cbot*ones(pgap,ngenes) , cbot*ones(pgap,ggap) , cbot*ones(pgap,size(cov_tmp,2)); ...
	    min(ctopl,max(cbotl,C(find(cov_ethn_b),:))) , cbot*ones(length(find(cov_ethn_b)),ggap) , cov_tmp(find(cov_ethn_b),:) ; ...
	     cbot*ones(pgap,ngenes) , cbot*ones(pgap,ggap) , cbot*ones(pgap,size(cov_tmp,2)); ...
	    min(ctopl,max(cbotl,C(find(cov_ethn_h),:))) , cbot*ones(length(find(cov_ethn_h)),ggap) , cov_tmp(find(cov_ethn_h),:) ; ...
	     cbot*ones(pgap,ngenes) , cbot*ones(pgap,ggap) , cbot*ones(pgap,size(cov_tmp,2)); ...
	    min(ctopl,max(cbotl,C(find(cov_ethn_c),:))) , cbot*ones(length(find(cov_ethn_c)),ggap) , cov_tmp(find(cov_ethn_c),:) ; ...
	 ],[cbot,ctop]); axis off;
title('sort by ethn (o,b,h,c)');
end;%plot_flag=0;

disp(sprintf('  '));
disp(sprintf(' While there are many different directions we can go with this data, '));
disp(sprintf(' for this example we ignore many of the possible covariates, '));
disp(sprintf(' electing to use only gender and ethnicity (caucasian vs non-caucasian). '));
cov_normal = cov_dss_no ;
cov_cancer = cov_dss_de ;
cov_ethn_x = cov_ethn_o | cov_ethn_b | cov_ethn_h ;
disp(sprintf(' We set up a covariate-matrix (cov_mat_use) which codes for these covariates. '));
disp(sprintf(' This cov_mat_use has 2 columns (one for each covariate category), '));
disp(sprintf(' each column listing the binary values associated with that covariate category. '));
disp(sprintf(' For this particular example each row of cov_mat_use takes the form of: '));
disp(sprintf(' [0,0] for female noncaucasian patients'));
disp(sprintf(' [1,0] for male noncaucasian patients'));
disp(sprintf(' [0,1] for female caucasian patients'));
disp(sprintf(' [1,1] for male caucasian patients'));
% setting up covariate matrix; 
cov_mat_use = [cov___male , cov_ethn_c];
disp(sprintf(' We also set up cov_cat_use which is a single column vector, '));
disp(sprintf(' each entry of which codes for the covariates of that patient. '));
disp(sprintf(' For this particular example each entry of cov_cat_use takes the form of: '));
disp(sprintf(' 1 + [0,0]*[1;2] = 1 for female noncaucasian patients'));
disp(sprintf(' 1 + [1,0]*[1;2] = 2 for male noncaucasian patients'));
disp(sprintf(' 1 + [0,1]*[1;2] = 3 for female caucasian patients'));
disp(sprintf(' 1 + [1,1]*[1;2] = 4 for male caucasian patients'));
cov_cat_use = (1 + 1*cov___male + 2*cov_ethn_c);
% to check cov_cat vs covariates ;
%figure;cla;[cv,ci] = sort(cov_cat_use);imagesc([cov_mat_use(ci,:) , cov_cat_use(ci)]);title('cov mat vs cov cat');

disp(sprintf('  '));
disp(sprintf(' Now we are in a position to normalize our original data (i.e., B) in preparation for biclustering '));
disp(sprintf(' The first kind of normalization we will implement is the simplest:  '));
disp(sprintf(' normalization of each gene independently across all patients (cases and controls). '));
disp(sprintf(' We then store all these (normalized) genes -- denoted as D --  for use in our biclustering algorithm. '));
disp(sprintf(' This particular file is called %s_n0x.mat, ',fpre_put));
disp(sprintf(' and also includes the number of patients (npats), number of genes (ngenes), '));
disp(sprintf(' patient-ids (Pnames), gene-ids (Gnames), the gene-indices to bicluster (g_ij), '));
disp(sprintf(' the case-patient-indices (d_ij), the control-patient-indices (x_ij), '));
disp(sprintf(' and the covariate data (cov_mat and cov_cat). '));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n0 normalization ; 
% we normalize each gene separately (considering all patients) ,
% then retain all genes ; 
% Note that the matrix "D" will be the MATRIX_ORIG for bfilter later on.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E = transpose((B)); D = zeros(size(E));
for ng = 1:ngenes;
if (mod(ng,1024)==0); disp(sprintf(' %% gene %d out of %d',ng,ngenes)); end;
tmp = E(find(cov_cancer | cov_normal),ng); tmp = (tmp-median(tmp))/std(tmp); D(find(cov_cancer | cov_normal),ng) = tmp;
end;%for ng = 1:ngenes;
% Here we retain all genes ;
g_ij_dxseg=[];
g_ij = setdiff(1:ngenes,g_ij_dxseg);
% storing indices that correspond to cases ;
d_ij = find(cov_cancer);
% and indices that correspond to controls. ;
x_ij = find(cov_normal);
cov_mat = cov_mat_use; cov_cat = cov_cat_use;
save(sprintf('%s_n0x.mat',fpre_put),'npats','ngenes','Pnames','Gnames','D','g_ij','d_ij','x_ij','cov_mat','cov_cat');

disp(sprintf('  '));
disp(sprintf(' The second kind of normalization we will implement again involves normalizing each gene independently. '));
disp(sprintf(' However, this time we throw out all the genes that are, by themselves, strongly differentially expressed '));
disp(sprintf(' across the cases, relative to the controls. '));
disp(sprintf(' Put another way, we choose to deliberately exclude the genes that are, say, highly expressed across the majority '));
disp(sprintf(' of cases and under-expressed across the majority of controls (or vice versa) '));
disp(sprintf(' In other words, we are deliberately ignoring the genes that most studies focus on! '));
disp(sprintf(' The main reason for this approach is so that we can search through the often ''discarded'' genes for biclusters. '));
disp(sprintf(' This particular file is called %s_n1x.mat, ',fpre_put));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n1 normalization ; 
% we normalize each gene separately (considering all patients) , 
% then retain genes which do not discriminate cases vs ctrls 
% (retaining those with p<0.15);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E = transpose((B)); D = zeros(size(E));
disp(sprintf('  '));
disp(sprintf(' First we perform the normalization: '));
for ng = 1:ngenes;
if (mod(ng,1024)==0); disp(sprintf(' %% processing gene %d out of %d',ng,ngenes)); end;
tmp = E(find(cov_cancer | cov_normal),ng); tmp = (tmp-median(tmp))/std(tmp); D(find(cov_cancer | cov_normal),ng) = tmp;
end;%for ng = 1:ngenes;
disp(sprintf('  '));
disp(sprintf(' Now we run through each of the genes, measuring whether or not each gene is differentially-expressed. '));
disp(sprintf(' We use the covariate-corrected auc of the case-values vs the control-values as a measure of how strongly '));
disp(sprintf(' any particular gene is differentially-expressed. '));
disp(sprintf(' Specifically, for each gene, we measure both the average-cauc and the minimum-cauc (across covariate categories). '));
disp(sprintf(' Histograms of these average-caucs (top left) and minimum-caucs (top right) are plotted in a figure. '));
cauc_avg_h0=zeros(ngenes,1); cauc_min_h0=zeros(ngenes,1);
for ng=1:ngenes;
if (mod(ng,1024)==0); disp(sprintf(' %% processing gene %d out of %d',ng,ngenes)); end;
[cauc_avg_h0(ng),cauc_min_h0(ng)] = get_corrected_auc(D(find(cov_cancer),ng),D(find(cov_normal),ng),cov_cat_use(find(cov_cancer),1),cov_cat_use(find(cov_normal,1)));
end;%for ng=1:ngenes;
subplot(2,2,1);hist(cauc_avg_h0,linspace(0,1,128)); xlabel('auc'); ylabel('# of genes'); title('histogram of cauc average'); subplot(2,2,2);hist(cauc_min_h0,linspace(0,1,128)); xlabel('auc'); ylabel('# of genes'); title('histogram of cauc minimum');
disp(sprintf('  '));
disp(sprintf(' To get a good idea of which genes are strongly differentially-expressed, we also calculate a distribution of cauc values '));
disp(sprintf(' under the ''covariate respecting'' label-shuffled null hypothesis. '));
disp(sprintf(' Specifically, this label-shuffled null hypothesis involves randomly interchanging case-control labels '));
disp(sprintf(' between patients of the same covariate type, '));
disp(sprintf(' For example, for each shuffled trial we randomly permute the case-control labels '));
disp(sprintf(' first within the female noncaucasian patients, '));
disp(sprintf(' then again within the male noncaucasian patents, '));
disp(sprintf(' then again within the female caucasian patients, '));
disp(sprintf(' and then finally within the male caucasian patients. '));
disp(sprintf(' We typically perform 32 permutations in total, although this number can certainly be increased. '));
%nprm = 32; 
clear nprm; nprm = input(' Input number of permutations? [default = 32]: '); if isempty(nprm); nprm=32; end;
disp(sprintf(' As we calculate the caucs for these permutations, we add to the figure above by '));
disp(sprintf(' showing the histograms of the average-cauc (lower left) and min-cauc (lower right) for each permutation. '));
cauc_avg_h1=zeros(ngenes,nprm); cauc_min_h1=zeros(ngenes,nprm);
for np=1:nprm;
for ns=max(1,min(cov_cat_use)):max(cov_cat_use);
%disp(sprintf(' %% permuting labels within cov_cat_use %d',ns));
std_d_use_ij{ns} = intersect(find(cov_cancer),find(cov_cat_use==ns));
std_x_use_ij{ns} = intersect(find(cov_normal),find(cov_cat_use==ns));
if (length(std_d_use_ij{ns})>0 & length(std_x_use_ij{ns})>0);
std_z_use_ij = [std_d_use_ij{ns} ; std_x_use_ij{ns}];
prm = randperm(length(std_z_use_ij));
std_d_use_ij{ns} = std_z_use_ij(prm(1:length(std_d_use_ij{ns})));
std_x_use_ij{ns} = std_z_use_ij(prm(length(std_d_use_ij{ns}) + (1:length(std_x_use_ij{ns}))));
end;%if (length(std_d_use_ij{ns})>0 & length(std_x_use_ij{ns})>0);
end;%for ns=max(1,min(cov_cat_use)):max(cov_cat_use);
ptmp_cancer = []; ptmp_cancer_cov_cat_use = [];
ptmp_normal = []; ptmp_normal_cov_cat_use = [];
for ns=max(1,min(cov_cat_use)):max(cov_cat_use);
tmp_l = length(std_d_use_ij{ns}); ptmp_cancer = [ptmp_cancer;std_d_use_ij{ns}]; ptmp_cancer_cov_cat_use = [ptmp_cancer_cov_cat_use;ns*ones(tmp_l,1)];
tmp_l = length(std_x_use_ij{ns}); ptmp_normal = [ptmp_normal;std_x_use_ij{ns}]; ptmp_normal_cov_cat_use = [ptmp_normal_cov_cat_use;ns*ones(tmp_l,1)];
end;%for ns=max(1,min(cov_cat_use)):max(cov_cat_use);
ptmp = randperm(npats);
for ng=1:ngenes;
if (mod(ng,1024)==0); disp(sprintf(' %% processing np %d/%d gene %d out of %d',np,nprm,ng,ngenes)); end;
[cauc_avg_h1(ng,np),cauc_min_h1(ng,np)] = get_corrected_auc(D(ptmp_cancer,ng),D(ptmp_normal,ng),ptmp_cancer_cov_cat_use,ptmp_normal_cov_cat_use);
end;%for ng=1:ngenes;
subplot(2,2,3);hist(cauc_avg_h1(:,np),linspace(0,1,128));
subplot(2,2,4);hist(cauc_min_h1(:,np),linspace(0,1,128));
drawnow;
end;%for np=1:nprm;
disp(sprintf('  '));
disp(sprintf(' Now we use this label-shuffled distribution to determine the 15th and 85th percentiles (across genes) for the '));
disp(sprintf(' measures of average-cauc and min-cauc. '));
disp(sprintf(' We use these percentiles as bounds to classify each gene in our original data as either: '));
disp(sprintf(' (a) possibly significantly differentially-expressed -- with either an average-cauc or min-cauc '));
disp(sprintf('     lying outside these percentile bounds, or'));
disp(sprintf(' (b) definitely not significantly differentially-expressed -- with an average-cauc and min-cauc '));
disp(sprintf('     both lying within the percentile bounds. '));
disp(sprintf(' The genes that fall into category-b are those that are usually discarded during a conventional analysis. '));
disp(sprintf(' In this situation we actually retain only the category-b genes, eliminating the genes in category-a. '));
pcutoff_avg = 0.15; plo_avg = prctile([cauc_avg_h1(:)],100*(pcutoff_avg)); phi_avg = prctile([cauc_avg_h1(:)],100*(1-pcutoff_avg));
pcutoff_min = 0.15; plo_min = prctile([cauc_min_h1(:)],100*(pcutoff_min)); phi_min = prctile([cauc_min_h1(:)],100*(1-pcutoff_min));
g_ij_dxseg_avg=find(cauc_avg_h0<plo_avg | cauc_avg_h0>phi_avg);
g_ij_dxseg_min=find(cauc_min_h0<plo_min | cauc_min_h0>phi_min);
g_ij_dxseg = union(g_ij_dxseg_avg,g_ij_dxseg_min);
% removing strongly case-control-dependent genes ;
g_ij = setdiff(1:ngenes,g_ij_dxseg);
disp(sprintf(' With the particular label-shuffled distribution generated above, we end up eliminating %d category-a genes',length(g_ij_dxseg)));
disp(sprintf(' and retaining %d category-b genes. ',length(g_ij)));
% storing indices that correspond to cases ;
d_ij = find(cov_cancer);
% and indices that correspond to controls. ;
x_ij = find(cov_normal);
cov_mat = cov_mat_use; cov_cat = cov_cat_use; 
save(sprintf('%s_n1x.mat',fpre_put),'npats','ngenes','Pnames','Gnames','D','g_ij','d_ij','x_ij','cov_mat','cov_cat');

disp(sprintf('  '));
disp(sprintf(' Finally, the third kind of normalization we will implement involves normalizing each gene across each ''group'' of patients, '));
disp(sprintf(' where the groups are determined by the case-control status and covariate categories. '));
disp(sprintf(' In this case we have three categories for each patient: case-control, male-female, and caucasian-noncaucasian; '));
disp(sprintf(' these three categories imply 2^3=8 groups of patients, each of which are normalized separately. '));
disp(sprintf(' This particular file is called %s_n2x.mat, ',fpre_put));
disp(sprintf(' Note that this particular normalization is designed to eliminate any differential expression that may exist '));
disp(sprintf(' across case-control status or across covariate categories. In other words, the only patterns that can remain '));
disp(sprintf(' after this normalization must involve co-expression within subsets of patients of a given group, relative to the mean of that group. '));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n2 normalization ;
% we normalize each gene across the different subgroups of patients, ;
% each subgroup associated with either case or ctrl ;
% as well as each unique combination of covariates ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

seg_cancer = find(cov_cancer);
seg_normal = find(cov_normal);
seg___male = find(cov___male);
seg_female = find(cov_female);
seg_ethn_c = find(cov_ethn_c);
seg_ethn_x = find(cov_ethn_x);
seg_cancer___male_ethn_c = find(cov_cancer & cov___male & cov_ethn_c);
seg_cancer_female_ethn_c = find(cov_cancer & cov_female & cov_ethn_c);
seg_cancer___male_ethn_x = find(cov_cancer & cov___male & cov_ethn_x);
seg_cancer_female_ethn_x = find(cov_cancer & cov_female & cov_ethn_x);
seg_normal___male_ethn_c = find(cov_normal & cov___male & cov_ethn_c);
seg_normal_female_ethn_c = find(cov_normal & cov_female & cov_ethn_c);
seg_normal___male_ethn_x = find(cov_normal & cov___male & cov_ethn_x);
seg_normal_female_ethn_x = find(cov_normal & cov_female & cov_ethn_x);
E = transpose((B)); D = zeros(size(E));
for ng = 1:ngenes;
if (mod(ng,1024)==0); disp(sprintf(' %% gene %d out of %d',ng,ngenes)); end;
tmp = E(seg_cancer___male_ethn_c,ng); tmp = (tmp-median(tmp))/std(tmp); D(seg_cancer___male_ethn_c,ng) = tmp;
tmp = E(seg_cancer_female_ethn_c,ng); tmp = (tmp-median(tmp))/std(tmp); D(seg_cancer_female_ethn_c,ng) = tmp;
tmp = E(seg_cancer___male_ethn_x,ng); tmp = (tmp-median(tmp))/std(tmp); D(seg_cancer___male_ethn_x,ng) = tmp;
tmp = E(seg_cancer_female_ethn_x,ng); tmp = (tmp-median(tmp))/std(tmp); D(seg_cancer_female_ethn_x,ng) = tmp;
tmp = E(seg_normal___male_ethn_c,ng); tmp = (tmp-median(tmp))/std(tmp); D(seg_normal___male_ethn_c,ng) = tmp;
tmp = E(seg_normal_female_ethn_c,ng); tmp = (tmp-median(tmp))/std(tmp); D(seg_normal_female_ethn_c,ng) = tmp;
tmp = E(seg_normal___male_ethn_x,ng); tmp = (tmp-median(tmp))/std(tmp); D(seg_normal___male_ethn_x,ng) = tmp;
tmp = E(seg_normal_female_ethn_x,ng); tmp = (tmp-median(tmp))/std(tmp); D(seg_normal_female_ethn_x,ng) = tmp;
end;%for ng = 1:ngenes;
disp(sprintf('  '));
disp(sprintf(' Now, by construction, there should be no differentially expressed genes. '));
disp(sprintf(' We can demonstrate this by measuring the average-cauc and min-cauc of the normalized data. '));
plot_flag = input(' Plot histogram of average-cauc and min-cauc? [default = 0]: '); if isempty(plot_flag); plot_flag=0; end;;  
if plot_flag;
cauc_avg_h0=zeros(ngenes,1); cauc_min_h0=zeros(ngenes,1);
for ng=1:ngenes;
if (mod(ng,1024)==0); disp(sprintf(' %% gene %d out of %d',ng,ngenes)); end;
[cauc_avg_h0(ng),cauc_min_h0(ng)] = get_corrected_auc(D(find(cov_cancer),ng),D(find(cov_normal),ng),cov_cat_use(find(cov_cancer),1),cov_cat_use(find(cov_normal),1));
end;%for ng=1:ngenes;
subplot(1,2,1);hist(cauc_avg_h0,linspace(0,1,128)); subplot(1,2,2);hist(cauc_min_h0,linspace(0,1,128)); drawnow();
end;%if plot_flag;
disp(sprintf(' Consequently, we remove no genes. '));
g_ij_dxseg = [];
g_ij = setdiff(1:ngenes,g_ij_dxseg);
% storing indices that correspond to cases ;
d_ij = find(cov_cancer);
% and indices that correspond to controls. ;
x_ij = find(cov_normal);
cov_mat = cov_mat_use; cov_cat = cov_cat_use; 
save(sprintf('%s_n2x.mat',fpre_put),'npats','ngenes','Pnames','Gnames','B','D','g_ij','d_ij','x_ij','cov_mat','cov_cat');
