function idx = num2idx(num_state, base_vector)
%IDX2NUM Summary of this function goes here
%   Detailed explanation goes here
    digit = 2;
    idx = num_state(1);
    while (digit <= size(base_vector,2) )
        idx = idx+ num_state(digit)*prod(base_vector(1,1:(digit-1)));
        digit = digit + 1;
    end
end

