function h = plot_cstl(constellation)
    M = length(constellation);
    h = plot(constellation, '.r', 'MarkerSize', 20);
    xlim([-1.5 1.5]);
    ylim([-1.5 1.5]);
    labels = cellstr(dec2bin(0:M-1));
    text(real(constellation),imag(constellation)+0.1,labels)
end