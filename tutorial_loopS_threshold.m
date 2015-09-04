function [output,ik_out,ij_out] = tutorial_loopS_threshold(LR2_ra,Lx2_ra,rthr,ar_min,ar_max,plot_flag);
% This function finds the corner-indices (ij_out,ik_out) of the array LR2_ra ;
% that satisfy the following constraints: ; 
% 1. the product ik_out*ij_out should be as large as possible (i.e., maximum area) ; 
% 2. the entry Lx2_ra(ij_out,ik_out) should be at a maximum (i.e., the full number of covariate-categories) ;
% 3. the aspect ratio (ij_out/nrows)/(ik_out/ncols) should lie in between ar_min and ar_max ;
% 4. the value LR2_ra(ij_out,ik_out) should be at least a threshold-value LR2_thr. ;
% Within the above constraints, the threshold-value LR2_thr is determined using the relative threshold rthr (an input). ;
% specifically, LR2_thr = LR2_min + rthr*(LR2_max - LR2_min) ;
% where LR2_min and LR2_max are the smallest and largest values of LR2_ra, respectively. ;
% The inputs to this function are:
% LR2_ra: an array of values (e.g., average row-scores) to threshold in constraint #4. ;
% Lx2_ra: an array of values (e.g., number of covariate-categories remaining) to use in constraint #2. ;
% rthr: a real number between 0 and 1 to use in constraint #4. ;
% ar_min, ar_max; positive real numbers to use in constraint #3. ;
% plot_flag: an integer (0 or 1); if this is set to 1 the output is plotted. ;

if (nargin<6); plot_flag=0; end;

[nrows,ncols] = size(LR2_ra);
output = zeros(1,1);
D = LR2_ra; 

[ij,ik] = find(Lx2_ra==max(Lx2_ra(:)));
LR2_max = 0;
for nl=1:length(ij);
LR2_max = max(LR2_max,LR2_ra(ij(nl),ik(nl)));
LR2_min = max(0,min(LR2_max,LR2_ra(ij(nl),ik(nl))));
LR2_thr = LR2_min + rthr*(LR2_max - LR2_min);
end;%for nl=1:length(ij);
[ij,ik] = find(LR2_ra>LR2_thr & Lx2_ra==max(Lx2_ra(:)));
loc_valid = find((ij/nrows)./(ik/ncols) > ar_min & (ij/nrows)./(ik/ncols) < ar_max);
if isempty(loc_valid); loc_valid=[1]; ij=1;ik=1; end;
output=0; if (length(loc_valid)>0); [output,loc] = max(ij(loc_valid).*ik(loc_valid)); end;%if (length(loc_valid)>0);

if plot_flag;
figure;cla;
if (length(loc_valid)>0);
hold on;
imagesc(LR2_ra(end:-1:1,:),[LR2_min,LR2_max]); colorbar;
for nl=1:length(loc_valid);
plot(ik(loc_valid(nl)),size(LR2_ra,1)-ij(loc_valid(nl)),'w.');
end;%for nl=1:length(loc_valid);
plot(ik(loc_valid(loc)),size(LR2_ra,1)-ij(loc_valid(loc)),'ko','Markersize',35);
plot(ik(loc_valid(loc)),size(LR2_ra,1)-ij(loc_valid(loc)),'kx','Markersize',35);
hold off;
xlim([1,size(LR2_ra,2)]); xlabel('columns n');
ylim([1,size(LR2_ra,1)]); ylabel('rows m');
set(gca,'Xtick',[ik(loc_valid(loc))],'Xticklabel',[ik(loc_valid(loc))],'Ytick',[size(LR2_ra,1)-ij(loc_valid(loc))],'Yticklabel',[ij(loc_valid(loc))]);
title(sprintf('valid areas in white; maximum area at ''x'': relative threshold %0.2f, abs threshold %0.2f',rthr,LR2_thr));
end;%if (length(loc_valid)>0);
end;%  if plot_flag;

ik_out = ik(loc_valid(loc));
ij_out = ij(loc_valid(loc));
%disp(sprintf(' %% ncov max %d, using rthr %0.3f = %0.2f [%0.2f,%0.2f] --> [%d-x-%d]=%d',max(Lx2_ra(:)),LR2_thr,rthr,LR2_min,LR2_max,ij_out,ik_out,output));
