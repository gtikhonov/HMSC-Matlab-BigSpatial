function [ A ] = insertRow( A, v, pos )
%INSERTROW Summary of this function goes here
%   Detailed explanation goes here
   if pos > 1
      if pos < size(A,1)+1 
         A = [A(1:(pos-1),:); v; A(pos:end,:)];
      else
         A = [A; v];
      end
   else
      A = [v; A];
   end
end

