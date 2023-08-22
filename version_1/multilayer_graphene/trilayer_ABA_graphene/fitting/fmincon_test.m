x0 = [-1,2];
p1 = 2; p2 =2;

%优化求解器选项配置
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
%新手可将 options  = []

%优化求解器返回最优解
x = fmincon(@fun,x0,[],[],[],[],[],[],@nonlcon,options,p1,p2); %返回优化后变量

function [f,gradf] = fun(x,p1,p2)
    f = p1*x(1)^2 + p2*x(2)^2;
    gradf = [2*p1*x(1);2*p2*x(2)];
end

% function [c,ceq] = nonlcon(x, p1, p2)
%     c = x(1)^2+2*x(1)*x(2)-5;
%     ceq = [];
% end

function [c,ceq] = nonlcon(x, p1, p2)
    c(1) = x(1)^2-x(2);
    c(2) = x(1)^2+2*x(1)*x(2)-5;
    ceq = x(1)-x(2)^2;
end

fun = @(x)100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
x0 = [0.5,0];
A = [1,2];
b = 1;
Aeq = [2,1];
beq = 1;
x = fmincon(fun,x0,A,b,Aeq,beq);