function output = unionall(input);
% unions a 1-dimensional cell_array of input;
l0=length(input); output = [];
for nl=1:l0; output = union(output,input{nl}); end;
