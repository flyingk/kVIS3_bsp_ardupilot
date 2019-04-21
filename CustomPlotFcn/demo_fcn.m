function y = demo_fcn(data, operator)

dt = mean(diff(operator));

y = fdfilt(data, 2, dt);