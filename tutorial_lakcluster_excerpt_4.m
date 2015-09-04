ncovs_rem = sum(cov_bother_flag_AAAA)+sum(cov_bother_flag_AZZA);
nrows_rem = nrows_rem - length(tmp_rij_drop_full);
ncols_rem = ncols_rem - length(tmp_cij_drop_full);
A = A(tmp_rij_keep_full,tmp_cij_keep_full);
Z = Z(:,tmp_cij_keep_full);
T = T(tmp_rij_keep_full,:);
S = S(:,:);
for ncov=0:ncovs-1;
%disp(sprintf(' %% ncov %d: lengths %d,%d,%d,%d',ncov,length(cov_A_up_{1+ncov}),length(cov_A_dn_{1+ncov}),length(cov_Z_up_{1+ncov}),length(cov_Z_dn_{1+ncov})));
cov_A_up_{1+ncov} = find(T(:,1+1+ncov)==1);
cov_A_dn_{1+ncov} = find(T(:,1+1+ncov)==0);
cov_Z_up_{1+ncov} = find(S(:,1+1+ncov)==1);
cov_Z_dn_{1+ncov} = find(S(:,1+1+ncov)==0);
%disp(sprintf(' %% ncov %d: lengths %d,%d,%d,%d',ncov,length(cov_A_up_{1+ncov}),length(cov_A_dn_{1+ncov}),length(cov_Z_up_{1+ncov}),length(cov_Z_dn_{1+ncov})));
end;%for ncov=0:ncovs-1;