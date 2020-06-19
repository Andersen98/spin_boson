function H_spin = H_s(bath_elements,epsilon)
%H_S Summary of this function goes here
%   Detailed explanation goes here
bath_space = prod(bath_elements(3,:));
bath_one = sparse(eye(bath_space));
mat1 = kron(sparse([0,1;1,0]),bath_one);
mat2 = kron(epsilon.*sparse([1,0;0,1]),bath_one);
H_spin = mat1+mat2;
end

