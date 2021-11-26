  function [mu,sigma] = eval_value_function(value_function,input)
% EVAL_VALUE_FUNCTION   Evaluation of value function.
%    [MU,SIGMA] = EVAL_VALUE_FUNCTION(VALUE_FUNCTION,INPUT,OUTPUT)
%    evaluates the value function at INPUT 
%    
%   As output, you will find:
%    MU the evalution at INPUT
%
%   See also FIT_VALUE_FUNCTION.



switch value_function.value_function_type
    case 'gp'
[mu,sigma]=value_function.gp_model(input,1);
    case 'nn'
%         if size(input,2)>size(value_function.nn.Y,2)
%          mu  = predict(value_function.nn,input(:,1));
%         else
%         mu  = predict(value_function.nn,input);
%         end
[mu,sigma]=value_function.nn_model(input);

        sigma = mu.*0;
end

  end