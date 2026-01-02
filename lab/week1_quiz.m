% impulse sequence
[x,n]=impseq(1,-2,5);
disp('impulse sequence');
disp([x;n]);

% unit step sequence
[x,n]=stepseq(1,-2,5);
disp('unit step sequence');
disp([x;n]);

% 

%impulse sequence function
function [x,n] = impseq(n0,n1,n2)
    n=[n1:n2];
    x=[(n-n0)==0];
end

% unit step sequence function
function [x,n] = stepseq(n0,n1,n2)
    n=[n1:n2];
    x=[(n-n0)>=0];
end
