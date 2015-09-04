ncovs_rem = sum(cov_bother_flag_AAAA)+sum(cov_bother_flag_AZZA);
nrows_rem = nrows_rem - length(tmp_rij_drop_full);
ncols_rem = ncols_rem - length(tmp_cij_drop_full);
if test_flag;
A = A(tmp_rij_keep_full,tmp_cij_keep_full);
Z = Z(:,tmp_cij_keep_full);
T = T(tmp_rij_keep_full,:);
S = S(:,:);
end;%if test_flag;
