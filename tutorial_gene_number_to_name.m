function [name_list,I__x,I_in,I_ls] = tutorial_gene_number_to_name(path_base,number_list,fname_out);
% converts number_list to name_list according to gene_entrez_symbol.txt ;

if ~(exist(sprintf('%s../gene_entrez_symbol.mat',path_base),'file'));
disp(sprintf('creating %s../gene_entrez_symbol.mat',path_base));
fid = fopen(sprintf('%s../gene_entrez_symbol.txt',path_base));
A = textscan(fid,'%s%s');
fclose(fid);
ngenes = length(A{1});
B = zeros(ngenes,1);
C = cell(ngenes,1);
for ng=1:ngenes;
B(ng) = str2num(A{1}{ng});
C{ng} = A{2}{ng};
end;%for ng=1:ngenes;
save(sprintf('%s../gene_entrez_symbol.mat',path_base),'ngenes','B','C');
end;%if ~(exist(sprintf('%s../gene_entrez_symbol.mat',path_base),'file'));
load(sprintf('%s../gene_entrez_symbol.mat',path_base),'ngenes','B','C');

[I__x,I_in,I_ls] = intersect(number_list,B,'stable'); name_list = C(I_ls);
if nargin>1;
fid = fopen(fname_out,'w');
for nl=1:length(I_in);
fprintf(fid,'%s\t%d\n',name_list{nl},number_list(I_in(nl)));
end;%for nl=1:length(I_in);
fclose(fid);
end;%if nargin>1;
