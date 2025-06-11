function FEMPiecewiseLinearErrConv()
% FEM_ERRORCONVERGENCE  Compute L2, energy, and max‑norm errors
% for the FEM approximation of u(x)=sin(pi x) on [0,1], with
% -u'' = π^2 sin(πx), u(0)=u(1)=0, using piecewise linear elements.
%
% We take h = 2^(-j), j=1..10, and plot each error vs h on a log‑log plot.

    clear; clc;

    Jmax = 10;
    H    = 2.^(- (1:Jmax));       % h values
    errL2  = zeros(1,Jmax);
    errE   = zeros(1,Jmax);
    errInf = zeros(1,Jmax);

    for j = 1:Jmax
        % h = H(j);
        x = (0:H(j):1)';          % nodes: length N+2
        N = length(x)-1;          % number of interior DOFs

        %--- assemble stiffness matrix A (NxN) ---%
        A_full = zeros(N+1,N+1);
        for i = 1:N
            h = x(i+1)-x(i);
            A_full(i,i) = A_full(i,i) + 1/h;
            A_full(i,i+1) = A_full(i,i+1) - 1/h;
            A_full(i+1,i) = A_full(i+1,i) - 1/h;
            A_full(i+1, i+1) = A_full(i+1, i+1) +  1/h;
        end
        % A(1,1) = A(1,1) + 0;
        % A(N+1,N+1) = A(N+1,N+1) + 0;
        A =  A_full(2:N, 2:N);
       % disp((A_full));
        

        %--- assemble load vector b (Nx1) by mass‑lumping ---%
        b_full = zeros(N+1,1);
        for i = 1:N
            h = x(i+1) - x(i);
            b_full(i)   = b_full(i)   + source(x(i))   * (h/2);
            b_full(i+1) = b_full(i+1) + source(x(i+1)) * (h/2);
        end
        % b(1) = b(1) + 0;
        % b(end) = b(end) + 0;
        b = b_full(2:end-1);

        %--- solve for interior values ---%
        u_int = A\b;
        u_h   = [0; u_int; 0];      % include boundary zeros

        %--- exact solution at nodes ---%
        u_ex_nodes = sin(pi*x);

        %--- (iii) max‑norm at nodes ---%
        errInf(j) = max(abs(u_ex_nodes - u_h));

        %--- (ii) energy norm: ||u-u_h||_E = sqrt((u_ex_int-u_h)'*A*(u_ex_int-u_h)) ---%
        w = u_ex_nodes(2:end-1) - u_int;
        errE(j) = sqrt( w' * A * w );

       %  --- (i) L2 norm via element‑wise Simpson’s rule ---%
        E2 = 0;
        for i = 1:length(x)-1
            xa = x(i); xb = x(i+1); xm = 0.5*(xa+xb);
            uh_a = u_h(i);    uh_b = u_h(i+1);    uh_m = 0.5*(uh_a+uh_b);
            ue_a = sin(pi*xa); ue_b = sin(pi*xb); ue_m = sin(pi*xm);
            % Simpson on [xa,xb] for (ue-uh)^2
            E2 = E2 + (xb-xa)/6 * ( (ue_a-uh_a)^2 ...
                                  + 4*(ue_m-uh_m)^2 ...
                                  + (ue_b-uh_b)^2 );
        end
        errL2(j) = sqrt(E2);
        % err_L2(j) = sqrt(h)*norm((u_ex_nodes(2:end-1)-u_int), 2);
    end

    %--- plot errors on a single log-log plot ---%
    figure; hold on; box on;
    loglog(H, errL2,  '-o','LineWidth',1.2);
    loglog(H, errE,   '-s','LineWidth',1.2);
    loglog(H, errInf,'-^','LineWidth',1.2);
    legend('||u-u_h||_{L^2}','||u-u_h||_E','max_i|u-u_h|','Location','SouthEast');
    xlabel('h'); ylabel('Error');
    title('Convergence of FEM for u(x)=\sin(\pi x)');
    grid on;

    %--- print observed rates ---%
    ratesL2  = log(errL2(1:end-1)./errL2(2:end)) ./ log(H(1:end-1)./H(2:end));
    ratesE   = log(errE(1:end-1)./errE(2:end))   ./ log(H(1:end-1)./H(2:end));
    ratesInf = log(errInf(1:end-1)./errInf(2:end)) ./ log(H(1:end-1)./H(2:end));
    fprintf(' h      rate_L2    rate_E    rate_max\n');
    for k=1:Jmax-1
      fprintf('%5.3f   %7.3f   %7.3f   %7.3f\n', H(k), ratesL2(k), ratesE(k), ratesInf(k));
    end
end

function y = source(x)
    % right‐hand side f(x) = pi^2 sin(pi x)
    y = pi^2 * sin(pi*x);
end
