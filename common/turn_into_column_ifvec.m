function y = turn_into_column_ifvec(x)

if any(size(x) == 1),
    y = x(:);
else
    y = x;
end