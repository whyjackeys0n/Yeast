clear;
clc;

syms x T R S P
PopPay = x^2*R+x*(1-x)*S+x*(1-x)*T+(1-x)^2*P;

PopPay_col = collect(PopPay,x);
PopPay_col_coef = coeffs(PopPay_col,x);
a = PopPay_col_coef(3);
b = PopPay_col_coef(2);
PopPay_max_x = collect(-b/2/a);

disp(['frequency = ',newline])
pretty(PopPay_max_x)