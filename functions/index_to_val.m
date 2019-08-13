function y = index_to_val(indices,vals)
% function y = index_to_val(indices,vals)
%converts index_nums by vals, ignoring nans
y = nan(size(indices));
inds = ~isnan(indices(:));
I = indices(inds);
y(inds) = vals(I);

