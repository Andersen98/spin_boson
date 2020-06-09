function coupling = coupling_generator(arg,params)
%Takes an argument to the type of coupling and params. 
%Spits out a list to describe coupling terms for each oscillator
if strcmp('none',arg)
    coupling = zeros(1,params.n_osc);
elseif strcmp('full',arg)
    coupling = ones(1,params.n_osc).*.1;
    coupling(1) = 1;
elseif strcmp('debye',arg)
    
else
    warning('error in coupling input')
    coupling = zeros(1,params.n_osc);
    
end

