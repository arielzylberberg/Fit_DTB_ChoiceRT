function y = symmetric_scale(xmax,dx)
% function y = symmetric_scale(xmax,dx)
% makes a set of values going from
% zero to the closest uppward rounded
% value to xmax. Then it adds negative 
% values

epsilon = 10^-10;
y = 0:dx:(xmax-epsilon+dx);
y = unique([-y,y]);
y = sort(y);