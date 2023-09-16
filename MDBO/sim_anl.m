function [x0, f0] = sim_anl(f, x0, lb, ub, Mmax, TolFun)
%%  模拟退火算法

% f ：适应度函数
% x0: 输入种群
% lb: 种群下边界
% ub: 种群上边界
% Mmax  : 最大温度
% TolFun: 优化变化容忍度
% x0: 输出优化后的种群
% f0: 输出优化后的种群的适应度值

% 计算适应度值
x = x0;
fx = feval(f, x);
f0 = fx;

% 模拟退火主要步骤
for m = 0: Mmax
    % 计算温度
    T = m / Mmax; 
    mu = 10^(T * 1000);

    for k=0:5
        dx = mu_inv(2 * rand - 1, mu) .* (ub - lb);
        x1 = x + dx;

        % 边界处理防止越界
        x1 = (x1 < lb) .* lb + (lb <= x1) .* (x1 <= ub) .* x1 + (ub < x1) .* ub;

        % 计算当前位置适应度值和适应度值偏差
        fx1 = feval(f, x1);
        df = fx1 - fx;

        % 如果df < 0则接受该解，如果大于0 则利用Metropolis准则进行判断是否接受       
        if (df < 0 || rand < exp(-T * df / (abs(fx) + eps) / TolFun)) == 1
            x = x1;
            fx = fx1;
        end

        % 判断当前解是否更优，更优则更新.       
        if (fx1 < f0) ==1
           x0 = x1;
           f0 = fx1;
        end   
    end
end
end