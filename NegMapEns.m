%% NegMapEns
% Computes the minimum negativity N of the normalized Choi operator of a
% channel, over all channels AB -> CD that map each state setrhoAB{i}
% into setsigmaCD{i}; if a value for epsilon is provided an approximate 
% mapping is considered.
% 
% required arguments:
% 
% setrhoAB:     Cell array containing the set of input states; all such
%               states should have the same input dimension dInTot
%
% dA:           Dimension of input system A
%
% setsigmaCD:   Cell array containing the set of output states; all such 
%               states should have the same output dimension dOutTot 
%
% dC:           Dimension of output system C
%
% optional arguments:
%
% epsilon:      When a value for epsilon is provided an approximate mapping
%               is considered, 
%                       such that: TraceNorm(W[rho] - sigma) <= epsilon
%
% Output:
%
% neg:          Minimum negativity of a normalised Choi state isomorphic to
%               a bipartite channel AB->CD that maps each input state
%               setrho{j} into the output state setsigma{j}
% W:            an optimal normalised Choi isomorphic state that achieves
%               the minimum negativity
%
% Additional package requirements:
%       + CVX     (http://cvxr.com/cvx/)
%
%       + Qetlab  (http://www.qetlab.com)
%             - Negativity
%             - Swap
%             - Partial Trace 
%             - Tensor
%             - TraceNorm
%
% authors: Marco Piani, Bintener Tom
%%-------------------------------------------------------------------------
function [neg,W_out] = NegMapEns(setrhoAB,dA,setsigmaCD,dC,varargin)

%% Calculate dimensions that are not given already

    %number of input states
    n = length(setrhoAB);
    %total input dimension (the dimension of each state setrhoAB{j}
    dInTot = length(setrhoAB{1});
    %total output dimension (the dimension of each state setsigmaAB{j}
    dOutTot = length(setsigmaCD{1});
    
%%--Error checking---------------------------------------------------------
     if n == length(setsigmaCD)
         for i = 2:n
             if length(setrhoAB{i}) ~= dInTot 
                 error('All input states should have the same dimension dInTot');
             elseif length(setsigmaCD{i}) ~= dOutTot
                 error('All output states should have the same dimension dOutTot');
             end 
         end
     else
         error('Number of input and output systems must be equal.')
     end
%%-------------------------------------------------------------------------     
     
    %dimension of input system B
    dB = dInTot/dA ;
    %dimension of output system D
    dD = dOutTot/dC ;
    % Set optional argument 
        if nargin == 4    
        elseif nargin == 5 
            if length(varargin{1}) == 1
                epsilon = varargin{1};               
            else
                error('Epsilon should be a number.')
            end
       else 
            error('Unexpected number of arguments.') 
       end


%% CVX part

cvx_begin sdp quiet
    %W is the Choi SDP variable
        variable W(dInTot*dOutTot,dInTot*dOutTot) hermitian semidefinite;

    %target funtion (Negativity) to be minimized
    %(We calculate the Negativity in the cut AC:BD so we have to reorder the
    %systems correctly using Swap)
        minimize Negativity(Swap(W,[2 3],[dC,dD,dA,dB]),dA*dC);
        
    % We impose necessary conditions
    subject to     

    %Trace-preservation of the channel is imposed via condition on 
    % partial trace of W; trace is taken over output systems C and D
        PartialTrace(W,[1,2],[dC,dD,dA,dB]) == eye(dInTot)/dInTot;
    
    %Imposing the mapping using W^(-1) of setrhoAB{i} into setsigmaCD{i}
            
        for i = 1:n   
            mapped = dInTot * PartialTrace(W*Tensor(eye(dOutTot),...
                     transpose(setrhoAB{i})),[3,4],[dC,dD,dA,dB]);
        
            if nargin == 5 && epsilon ~= 0
                TraceNorm(mapped-setsigmaCD{i}) <= epsilon; 
            else
                mapped == setsigmaCD{i};
           end
        end   
           
    %End of CVX    
    cvx_end

%% Determine the output

	%minimum negativity required
        neg = cvx_optval;
    %Choi of an optimal map
        W_out = W;
end

