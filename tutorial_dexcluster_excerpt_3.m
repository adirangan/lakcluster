%rdrop = r_rem{iteration}(tmp_rij(1:ceil(gamma_tmp_row*end))); cdrop = c_rem{iteration}(tmp_cij(1:ceil(gamma_tmp_col*end)));
%rkeep = r_rem{iteration}(tmp_rij(ceil(gamma_tmp_row*end)+1:end)); ckeep = c_rem{iteration}(tmp_cij(ceil(gamma_tmp_col*end)+1:end));
rdrop = tmp_rij_drop_full; cdrop = tmp_cij_drop_full;
% for each dropped row we modify X_T_all, and for each ncov we modify both X_T_xx_up(nc,ncov) and X_T_xx_dn(nc,ncov), and some of the X_T_up_up, X_T_up_dn, X_T_dn_up, X_T_dn_dn ;
for nr=1:length(rdrop);
X_T_all(:) = X_T_all(:) - transpose(X_P_all(rdrop(nr),:)); X_P_all(rdrop(nr),:)=0;
for ncov=1:ncovs;
X_T_xx_up(:,ncov) = X_T_xx_up(:,ncov) - transpose(X_P_xx_up{ncov}(rdrop(nr),:)); X_P_xx_up{ncov}(rdrop(nr),:)=0;
X_T_xx_dn(:,ncov) = X_T_xx_dn(:,ncov) - transpose(X_P_xx_dn{ncov}(rdrop(nr),:)); X_P_xx_dn{ncov}(rdrop(nr),:)=0;
tmp_ij = find(cov_A_up_orig_{ncov}==rdrop(nr));
if ~isempty(tmp_ij); 
X_T_up_up(:,ncov) = X_T_up_up(:,ncov) - transpose(X_P_up_up{ncov}(tmp_ij,:)); X_P_up_up{ncov}(tmp_ij,:)=0;
X_T_up_dn(:,ncov) = X_T_up_dn(:,ncov) - transpose(X_P_up_dn{ncov}(tmp_ij,:)); X_P_up_dn{ncov}(tmp_ij,:)=0;
cov_A_up_{ncov} = setdiff(cov_A_up_{ncov},rdrop(nr));
end;% if found;
tmp_ij = find(cov_A_dn_orig_{ncov}==rdrop(nr));
if ~isempty(tmp_ij); 
X_T_dn_up(:,ncov) = X_T_dn_up(:,ncov) - transpose(X_P_dn_up{ncov}(tmp_ij,:)); X_P_dn_up{ncov}(tmp_ij,:)=0;
X_T_dn_dn(:,ncov) = X_T_dn_dn(:,ncov) - transpose(X_P_dn_dn{ncov}(tmp_ij,:)); X_P_dn_dn{ncov}(tmp_ij,:)=0;
cov_A_dn_{ncov} = setdiff(cov_A_dn_{ncov},rdrop(nr));
end;% if found;
end;%for ncov=1:ncovs;
end;%for nr=1:length(rdrop);
% for each dropped col we modify X_P_all, and for each ncov we modify X_P_xx_up, X_P_xx_dn, and some of the X_P_up_up, X_P_up_dn, X_P_dn_up, X_P_dn_dn ;
for nc=1:length(cdrop);
%disp(sprintf('setting cdrop(%d)=%d to zero',nc,cdrop(nc)));
X_P_all(:,cdrop(nc)) = 0;
for ncov=1:ncovs;
X_P_xx_up{ncov}(:,cdrop(nc)) = 0;
X_P_xx_dn{ncov}(:,cdrop(nc)) = 0;
X_P_up_up{ncov}(:,cdrop(nc)) = 0;
X_P_up_dn{ncov}(:,cdrop(nc)) = 0;
X_P_dn_up{ncov}(:,cdrop(nc)) = 0;
X_P_dn_dn{ncov}(:,cdrop(nc)) = 0;
end;%for ncov=1:ncovs;
end;%for nc=1:length(cdrop);
