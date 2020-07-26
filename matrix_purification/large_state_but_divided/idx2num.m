function num_state = idx2num(idx, base_vector)
%IDX2NUM Summary of this function goes here
%   Detailed explanation goes here
    count = 1;
    num_state = zeros(1,size(base_vector,2));
    while (idx ~=0 )
        num_state(1,count) = mod(idx,base_vector(1,count));
        idx = floor(idx/base_vector(1,count));
        count = count + 1;
    end
end

