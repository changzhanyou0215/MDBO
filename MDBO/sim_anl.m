function [x0, f0] = sim_anl(f, x0, lb, ub, Mmax, TolFun)
%%  ģ���˻��㷨

% f ����Ӧ�Ⱥ���
% x0: ������Ⱥ
% lb: ��Ⱥ�±߽�
% ub: ��Ⱥ�ϱ߽�
% Mmax  : ����¶�
% TolFun: �Ż��仯���̶�
% x0: ����Ż������Ⱥ
% f0: ����Ż������Ⱥ����Ӧ��ֵ

% ������Ӧ��ֵ
x = x0;
fx = feval(f, x);
f0 = fx;

% ģ���˻���Ҫ����
for m = 0: Mmax
    % �����¶�
    T = m / Mmax; 
    mu = 10^(T * 1000);

    for k=0:5
        dx = mu_inv(2 * rand - 1, mu) .* (ub - lb);
        x1 = x + dx;

        % �߽紦����ֹԽ��
        x1 = (x1 < lb) .* lb + (lb <= x1) .* (x1 <= ub) .* x1 + (ub < x1) .* ub;

        % ���㵱ǰλ����Ӧ��ֵ����Ӧ��ֵƫ��
        fx1 = feval(f, x1);
        df = fx1 - fx;

        % ���df < 0����ܸý⣬�������0 ������Metropolis׼������ж��Ƿ����       
        if (df < 0 || rand < exp(-T * df / (abs(fx) + eps) / TolFun)) == 1
            x = x1;
            fx = fx1;
        end

        % �жϵ�ǰ���Ƿ���ţ����������.       
        if (fx1 < f0) ==1
           x0 = x1;
           f0 = fx1;
        end   
    end
end
end