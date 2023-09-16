function x = mu_inv(y, mu)

% 模拟退火产生新位置偏差
x = (((1 + mu).^abs(y) - 1) / mu) .* sign(y);

end