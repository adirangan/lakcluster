function output = intersectall(input);
% intersects a 1-dimensional cell_array of input;
l0=length(input); output = [];
for nl=1:l0; if (nl==1); output = input{nl}; else output = intersect(output,input{nl}); end; end;
