function y = positive_constrained(x)
    if x < 0
        y = exp(x);
    else
        y = x + 1;
    end
end

