%% prod_state 
% Computes the Kronecker product of two states |n> and |m> of dimension d.
%
% Required arguments:
%
% n:    Index of the first state in the standard basis.
% m:    Index of the first state in the standard basis.
%
% Optional argument
%
% d:    Dimension of the states, default is set to 2. 
%
% Output:
% 
% prod_state: The product state in the standard basis.
%
% authors: Bintener Tom
%%-------------------------------------------------------------------------
function[prod_state] = prod_state(n,m,varargin)
    % Determine the value of d and check that it is at least 2.
    switch nargin
        case 2 
            d = 2;
        case 3 
            if varargin{1} < 2
                error('d needs to be at least 2.')
            else
                d = varargin{1};
            end
        otherwise
            error('Unexpected number of arguments.')
    end
    
    % Verify that both n and m are strictly less than d.
     if n >= d
         error('n needs to be strictly less than d.')
     elseif m >= d
         error('m needs to be strictly less than d.')
     elseif n < 0 || m < 0
         error('n and m need to be greater or equal to zero.')
     else
        I = eye(d);
        prod_state = kron(I(:,n+1),I(:,m+1));
     end
end