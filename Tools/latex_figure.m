function [] = latex_figure()

    reset(groot)
    set(groot,'defaulttextinterpreter','latex')
    set(groot,'defaultaxesticklabelinterpreter','latex')
    set(groot,'defaultlegendinterpreter','latex')
    set(groot,'defaultAxesXgrid','on','defaultAxesYgrid',...
        'on','defaultAxesZgrid','on')
    set(groot,'defaultAxesGridLineStyle','-.')
end

