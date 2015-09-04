function [XR_AAAA_full,XC_AAAA_full,XR_AZZA_full,XC_AZZA_full,X__AnAt_full,X__AtAn_full,X__AnZt_full,X__ZtZn_full] = tutorial_lakcluster_update_0(ver_flag,tmp_rij_keep,tmp_rij_drop,tmp_cij_keep,tmp_cij_drop,A,Z,XR_AAAA_full,XC_AAAA_full,XR_AZZA_full,XC_AZZA_full,X__AnAt_full,X__AtAn_full,X__AnZt_full,X__ZtZn_full);
% updates row and column scores ;
% test with: ;
%{

  % make sure error==0;
  clear all;
  MA = 81;  NA = 140;
  MZ = 79;  NZ = NA;
  A = randn(MA,NA);  Z = randn(MZ,NZ);
  tmp_rij_drop = [5 2 7]; tmp_rij_keep = setdiff(1:MA,tmp_rij_drop);
  tmp_cij_drop = [9 17 11 23]; tmp_cij_keep = setdiff(1:NA,tmp_cij_drop);
  X__AnAt_full = A*transpose(A);
  X__AtAn_full = transpose(A)*A;
  X__AnZt_full = A*transpose(Z);
  X__ZtZn_full = transpose(Z)*Z;
  XR_AAAA_full = diag(X__AnAt_full*transpose(X__AnAt_full));
  XC_AAAA_full = diag(transpose(A)*X__AnAt_full*A);
  XR_AZZA_full = diag(X__AnZt_full*transpose(X__AnZt_full));
  XC_AZZA_full = diag(transpose(A)*X__AnZt_full*Z);
  ver_flag=1;
  if ver_flag==0;
  [XR_AAAA_full,XC_AAAA_full,XR_AZZA_full,XC_AZZA_full,X__AnAt_full,X__AtAn_full,X__AnZt_full,X__ZtZn_full] = tutorial_lakcluster_update_0(ver_flag,tmp_rij_keep,tmp_rij_drop,tmp_cij_keep,tmp_cij_drop,A,Z,XR_AAAA_full,XC_AAAA_full,XR_AZZA_full,XC_AZZA_full,X__AnAt_full,X__AtAn_full,X__AnZt_full,X__ZtZn_full);
  elseif ver_flag==1;
  [XR_AAAA_full,XC_AAAA_full,XR_AZZA_full,XC_AZZA_full] = tutorial_lakcluster_update_0(ver_flag,tmp_rij_keep,tmp_rij_drop,tmp_cij_keep,tmp_cij_drop,A,Z,XR_AAAA_full,XC_AAAA_full,XR_AZZA_full,XC_AZZA_full);
  end;% ver_flag;
  A_out = A(tmp_rij_keep,tmp_cij_keep); Z_out = Z(:,tmp_cij_keep);
  X__AnAt__out = A_out*transpose(A_out);
  X__AtAn__out = transpose(A_out)*A_out;
  X__AnZt__out = A_out*transpose(Z_out);
  X__ZtZn__out = transpose(Z_out)*Z_out;
  XR_AAAA__out = diag(X__AnAt__out*transpose(X__AnAt__out));
  XC_AAAA__out = diag(transpose(A_out)*X__AnAt__out*A_out);
  XR_AZZA__out = diag(X__AnZt__out*transpose(X__AnZt__out));
  XC_AZZA__out = diag(transpose(A_out)*X__AnZt__out*Z_out);
  if ver_flag==0;
  n__AnAt = log(norm(X__AnAt_full - X__AnAt__out)); n__AtAn = log(norm(X__AtAn_full - X__AtAn__out)); n__AnZt = log(norm(X__AnZt_full - X__AnZt__out)); n__ZtZn = log(norm(X__ZtZn_full - X__ZtZn__out)); 
  nR_AAAA = log(norm(XR_AAAA_full - XR_AAAA__out)); nC_AAAA = log(norm(XC_AAAA_full - XC_AAAA__out)); nR_AZZA = log(norm(XR_AZZA_full - XR_AZZA__out)); nC_AZZA = log(norm(XC_AZZA_full - XC_AZZA__out));
  disp(sprintf(' %% errors: exp([%.2f,%.2f,%.2f,%.2f ; %.2f,%.2f,%.2f,%.2f])',n__AnAt,n__AtAn,n__AnZt,n__ZtZn,nR_AAAA,nC_AAAA,nR_AZZA,nC_AZZA));
  elseif ver_flag==1;
  nR_AAAA = log(norm(XR_AAAA_full - XR_AAAA__out)); nC_AAAA = log(norm(XC_AAAA_full - XC_AAAA__out)); nR_AZZA = log(norm(XR_AZZA_full - XR_AZZA__out)); nC_AZZA = log(norm(XC_AZZA_full - XC_AZZA__out));
  disp(sprintf(' %% errors: exp([%.2f,%.2f,%.2f,%.2f])',nR_AAAA,nC_AAAA,nR_AZZA,nC_AZZA));
  end;% ver_flag;

  % test timing ;
  clear all;
  profile on;
  disp(sprintf(' %% initializing... '));
  MA = 32;  NA = 3100;
  MZ = 33;  NZ = NA;
  A = randn(MA,NA);  Z = randn(MZ,NZ);
  tmp_rij_drop = [5 2 7]; tmp_rij_keep = setdiff(1:MA,tmp_rij_drop);
  tmp_cij_drop = [9 17 11 23]; tmp_cij_keep = setdiff(1:NA,tmp_cij_drop);
  tic;
  X__AnAt_full = A*transpose(A);
  X__AtAn_full = transpose(A)*A;
  X__AnZt_full = A*transpose(Z);
  X__ZtZn_full = transpose(Z)*Z;
  XR_AAAA_full = diag(X__AnAt_full*transpose(X__AnAt_full));
  XC_AAAA_full = diag(transpose(A)*X__AnAt_full*A);
  XR_AZZA_full = diag(X__AnZt_full*transpose(X__AnZt_full));
  XC_AZZA_full = diag(transpose(A)*X__AnZt_full*Z);
  tini = toc;
  disp(sprintf(' %% init time %f',tini));
  max_iteration = 1;
  disp(sprintf(' %% running %d iterations',max_iteration));
  for (ver_flag=[1]);
  tic;
  for iteration=1:max_iteration;
  if ver_flag==0;
  [XR_AAAA__tmp,XC_AAAA__tmp,XR_AZZA__tmp,XC_AZZA__tmp,X__AnAt__tmp,X__AtAn__tmp,X__AnZt__tmp,X__ZtZn__tmp] = tutorial_lakcluster_update_0(ver_flag,tmp_rij_keep,tmp_rij_drop,tmp_cij_keep,tmp_cij_drop,A,Z,XR_AAAA_full,XC_AAAA_full,XR_AZZA_full,XC_AZZA_full,X__AnAt_full,X__AtAn_full,X__AnZt_full,X__ZtZn_full);
  elseif ver_flag==1;
  [XR_AAAA__tmp,XC_AAAA__tmp,XR_AZZA__tmp,XC_AZZA__tmp] = tutorial_lakcluster_update_0(ver_flag,tmp_rij_keep,tmp_rij_drop,tmp_cij_keep,tmp_cij_drop,A,Z,XR_AAAA_full,XC_AAAA_full,XR_AZZA_full,XC_AZZA_full);
  end;% ver_flag;
  end;% for iteration=1:max_iteration;
  tsum(1+ver_flag) = toc;
  end;% for (ver_flag=[0,1]);
  disp(sprintf(' %% tsum %f %f  tper %f %f ',tsum,tsum/max_iteration));  
  profile off;
  profile viewer;

%}

tmp_An = A(tmp_rij_keep,tmp_cij_keep);
tmp_Cn = A(tmp_rij_keep,tmp_cij_drop);
tmp_Bt = A(tmp_rij_drop,tmp_cij_keep);
tmp_Dn = A(tmp_rij_drop,tmp_cij_drop);
tmp_Zn = Z(:,tmp_cij_keep);
tmp_En = Z(:,tmp_cij_drop);
if ver_flag==0;
tmp_CnCt = tmp_Cn*transpose(tmp_Cn);
tmp_AnBn = tmp_An*transpose(tmp_Bt);
tmp_CnDt = tmp_Cn*transpose(tmp_Dn);
tmp_ABCD = tmp_AnBn + tmp_CnDt;
tmp_BnBt = transpose(tmp_Bt)*tmp_Bt;
tmp_AtCn = transpose(tmp_An)*tmp_Cn;
tmp_BnDn = transpose(tmp_Bt)*tmp_Dn;
tmp_ACBD = tmp_AtCn + tmp_BnDn;
X__AnAt_full = X__AnAt_full(tmp_rij_keep,tmp_rij_keep) - tmp_CnCt ;
XR_AAAA_full = XR_AAAA_full(tmp_rij_keep) - diag( X__AnAt_full*tmp_CnCt + tmp_CnCt*X__AnAt_full + tmp_CnCt*tmp_CnCt + tmp_ABCD*transpose(tmp_ABCD) ) ;
X__AtAn_full = X__AtAn_full(tmp_cij_keep,tmp_cij_keep) - tmp_BnBt ;
XC_AAAA_full = XC_AAAA_full(tmp_cij_keep) - diag( X__AtAn_full*tmp_BnBt + tmp_BnBt*X__AtAn_full + tmp_BnBt*tmp_BnBt + tmp_ACBD*transpose(tmp_ACBD) ) ;
tmp_CnEt = tmp_Cn*transpose(tmp_En);
X__AnZt_full = X__AnZt_full(tmp_rij_keep,:) - tmp_CnEt ;
XR_AZZA_full = XR_AZZA_full(tmp_rij_keep) - diag( X__AnZt_full*transpose(tmp_CnEt) + tmp_CnEt*transpose(X__AnZt_full) + tmp_CnEt*transpose(tmp_CnEt) ) ;
tmp_ZtEn = transpose(tmp_Zn)*tmp_En;
X__ZtZn_full = X__ZtZn_full(tmp_cij_keep,tmp_cij_keep) ;
XC_AZZA_full = XC_AZZA_full(tmp_cij_keep) - diag( tmp_BnBt*X__ZtZn_full + tmp_ACBD*transpose(tmp_ZtEn) ) ;
elseif ver_flag==1;
tmp_CtCn = transpose(tmp_Cn)*tmp_Cn;
tmp_AnBn = tmp_An*transpose(tmp_Bt);
tmp_CnDt = tmp_Cn*transpose(tmp_Dn);
tmp_ABCD = tmp_AnBn + tmp_CnDt;
tmp_BtBn = tmp_Bt*transpose(tmp_Bt);
tmp_AtCn = transpose(tmp_An)*tmp_Cn;
tmp_BnDn = transpose(tmp_Bt)*tmp_Dn;
tmp_ACBD = tmp_AtCn + tmp_BnDn;
tmp_CnEt = tmp_Cn*transpose(tmp_En);
tmp_ZtEn = transpose(tmp_Zn)*tmp_En;
tmp_EtEn = transpose(tmp_En)*tmp_En;
tmp_ZnBn = tmp_Zn*transpose(tmp_Bt);
XR_AAAA__tmp = zeros(length(tmp_rij_keep),1);
XR_AZZA__tmp = zeros(length(tmp_rij_keep),1);
for nl=1:length(tmp_rij_keep);
tmp_an = tmp_An(nl,:); tmp_at = transpose(tmp_an); tmp_cn = tmp_Cn(nl,:); tmp_ct = transpose(tmp_cn); 
XR_AAAA__tmp(nl) = (tmp_an*tmp_AtCn)*tmp_ct + tmp_cn*transpose(tmp_AtCn)*tmp_at + tmp_cn*tmp_CtCn*tmp_ct + tmp_ABCD(nl,:)*transpose(tmp_ABCD(nl,:));
XR_AZZA__tmp(nl) = (tmp_an*tmp_ZtEn)*tmp_ct + tmp_cn*transpose(tmp_ZtEn)*tmp_at + tmp_cn*tmp_EtEn*tmp_ct;
end;%for nl=1:length(tmp_rij_keep);
XR_AAAA_full = XR_AAAA_full(tmp_rij_keep) - XR_AAAA__tmp;
XR_AZZA_full = XR_AZZA_full(tmp_rij_keep) - XR_AZZA__tmp;
XC_AAAA__tmp = zeros(length(tmp_cij_keep),1);
XC_AZZA__tmp = zeros(length(tmp_cij_keep),1);
for nl=1:length(tmp_cij_keep);
tmp_an = tmp_An(:,nl); tmp_at = transpose(tmp_an); tmp_bt = tmp_Bt(:,nl); tmp_bn = transpose(tmp_bt); tmp_zn = tmp_Zn(:,nl); tmp_zt = transpose(tmp_zn);
XC_AAAA__tmp(nl) = (tmp_at*tmp_AnBn)*tmp_bt + (tmp_bn*transpose(tmp_AnBn))*tmp_an + (tmp_bn*tmp_BtBn)*tmp_bt + tmp_ACBD(nl,:)*transpose(tmp_ACBD(nl,:));
XC_AZZA__tmp(nl) = (tmp_bn*transpose(tmp_ZnBn))*tmp_zn + tmp_ACBD(nl,:)*transpose(tmp_ZtEn(nl,:));
end;%for nl=1:length(tmp_cij_keep);
XC_AAAA_full = XC_AAAA_full(tmp_cij_keep) - XC_AAAA__tmp;
XC_AZZA_full = XC_AZZA_full(tmp_cij_keep) - XC_AZZA__tmp;
 else; disp(sprintf(' %% Warning! wrong version %d within tutorial_lakcluster_update_1.m',ver_flag)); return;
end;%if ver_flag==1;
