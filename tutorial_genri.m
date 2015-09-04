function tutorial_genri(path_base,gen_fname,frwd_vs_back,fname_out);
% function tutorial_genri();
% designed to evalute gene-enrichments of gene-sets stored in files such as GSE3307_n2x_genes_frwd_1_glist.txt ;
% note that this uses seek.GeneEnrichTest. ;
% Must run something like: ;
% bash; cd /data/rangan/dir_bcc/dir_code_012015/hierarchical/; javac seek/GeneEnrichmentTest.java;
% before seek will function ;

bitj=16; 
loopS_flag=1;
prefix_base = sprintf('%s',gen_fname); prefix=prefix_base; 

if exist(fname_out,'file'); fid_out = fopen(fname_out,'a');
elseif ~(exist(fname_out,'file'));
 fid_out = fopen(fname_out,'w');
fprintf(fid_out,'bash ;\n');
fprintf(fid_out,'cd /data/rangan/dir_bcc/dir_code_021015/hierarchical/ ;\n');
fprintf(fid_out,'source init.sh ;\n');
fprintf(fid_out,'javac seek/GeneEnrichmentTest.java ;\n');
fprintf(fid_out,'java seek.GeneEnrichmentTest -h ;\n');
end;% if exist ;

NRUN=16;
for nrun=0:NRUN;
path_use = path_base; prefix = prefix_base;
if (nrun==0); path_plus_prefix = sprintf('%s%s_%s',path_use,prefix,frwd_vs_back);
 else path_plus_prefix = sprintf('%s%s_%s_%d',path_use,prefix,frwd_vs_back,nrun); end;%if (nrun==0);
glist_fname = sprintf('%s_glist.txt',path_plus_prefix);
if exist(glist_fname,'file');
disp(sprintf(' %% found file %s',glist_fname));
glist_tmp = textread(glist_fname);
genri_fname = sprintf('%s_genri.txt',path_plus_prefix);
command_string = sprintf('java seek.GeneEnrichmentTest -l /data/rangan/dir_bcc/dir_code_021015/hierarchical/enrichment/ -m 2 %s -s original -t %d > %s ;',glist_fname,min(2048,length(glist_tmp)),genri_fname);
disp(command_string); fprintf(fid_out,sprintf('%s\n',command_string)); 
command_string = sprintf('java seek.GeneEnrichmentTest -l /data/rangan/dir_bcc/dir_code_021015/hierarchical/enrichment/ -m 2 %s -s original -t %d >! %s ;',glist_fname,min(2048,length(glist_tmp)),genri_fname);
system(command_string);
genri_fname = sprintf('%s_gslim.txt',path_plus_prefix);
command_string = sprintf('java seek.GeneEnrichmentTest -l /data/rangan/dir_bcc/dir_code_021015/hierarchical/enrichment/ -m 2 %s -s experimental_bp_slim -t %d > %s ;',glist_fname,min(2048,length(glist_tmp)),genri_fname);
disp(command_string); fprintf(fid_out,sprintf('%s\n',command_string)); 
command_string = sprintf('java seek.GeneEnrichmentTest -l /data/rangan/dir_bcc/dir_code_021015/hierarchical/enrichment/ -m 2 %s -s experimental_bp_slim -t %d >! %s ;',glist_fname,min(2048,length(glist_tmp)),genri_fname);
system(command_string);
 else if (nrun==0); disp(sprintf(' %% cannot find file %s',glist_fname)); end;
end;%if exist(glist_fname,'file');
end;%for nrun=0:NRUN;

fclose(fid_out);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

