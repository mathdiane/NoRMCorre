function W = construct_weights(regular_grid,extended_grid)
%Inputs: extend_grid comes from the regular_grid extended by overlapping

%trapezoidal fx

%linspace(x1,x2,n): generates n points btw [x1,x2];
%xx = [a value close to 0 -> 1 (overlapping), many 1's for regular_grid
%size, 1 -> a value close to 0 (overlapping)]
%similarly for yy and zz
xx = [linspace(1/(regular_grid(1)-extended_grid(1)+1),1,regular_grid(1)-extended_grid(1)),ones(1,regular_grid(2)-regular_grid(1)+1),fliplr(linspace(1/(extended_grid(2)-regular_grid(2)+1),1,extended_grid(2)-regular_grid(2)))];
yy = [linspace(1/(regular_grid(3)-extended_grid(3)+1),1,regular_grid(3)-extended_grid(3)),ones(1,regular_grid(4)-regular_grid(3)+1),fliplr(linspace(1/(extended_grid(4)-regular_grid(4)+1),1,extended_grid(4)-regular_grid(4)))];
zz = [linspace(1/(regular_grid(5)-extended_grid(5)+1),1,regular_grid(5)-extended_grid(5)),ones(1,regular_grid(6)-regular_grid(5)+1),fliplr(linspace(1/(extended_grid(6)-regular_grid(6)+1),1,extended_grid(6)-regular_grid(6)))];

[XX,YY,ZZ] = meshgrid(xx,yy,zz);
% The grid vectors xx,yy,zz form the rows of XX, columns of YY, and pages of Z 
% respectively. XX is xx repeatedly stacked vertically (along 2-dim) for
% length(yy) times; YY is yy repeatedly stacked horizontally (along 1st dim) for
% length(xx) times; if zz==1; ZZ an length(xx) x length(yy) array that's
% full of 1's.
W = XX.*YY.*ZZ;