function myplot(X, Y, ylabel_s, title_str, n)
%PLOTWAV 此处显示有关此函数的摘要
%   此处显示详细说明

figure(n);
plot(X,Y);
axis xy
axis tight
ylabel(ylabel_s);
xlabel('Time (secs)');
title(title_str);
end

