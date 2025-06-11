# Continuous-Piecewise-Linear-FEM
A MATLAB implementation to compute and visualize the convergence of the finite element method (FEM) with piecewise linear elements for the dirichlet boundary‑value problem: -u'' = pi sin(pi x) in (0,1)

# Description 
This repository contains a MATLAB code that computes the L²-norm, energy norm, and maximum norm errors for the FEM approximation. It then plots these errors against mesh size  on a log–log plot and prints the observed convergence rates.

# Code Overview 
The main script FEMPiecewiseLinearErrConv performs the following steps: Loop over mesh refinements h = 2^{-j}, j= 1,..,J_max.
- Assemble the stiffness matrix via element contributions.
- Assemble the load vector by mass-lumping the source term f = sin(pi x).
- Solve the linear system for interior nodal values, enforcing zero Dirichlet boundary conditions.

Compute error norms:
 - L² norm via Simpson’s rule on each element.
 - Energy norm via Trapezoidal rule on each element. 
 - Max norm (∞-norm) at nodes.

Plot all three errors vs.  on a single log–log plot. Print the observed convergence rates for each norm.

# Main Function 
```matlab
function FEMPiecewiseLinearErrConv()
    clear; clc;

    Jmax   = 10;
    H      = 2.^(- (1:Jmax));       % mesh sizes h
    errL2  = zeros(1,Jmax);
    errE   = zeros(1,Jmax);
    errInf = zeros(1,Jmax);

    for j = 1:Jmax
        x      = (0:H(j):1)';          % global nodes
        N      = length(x)-1;          % # intervals

        % Assemble stiffness matrix A
        A_full = zeros(N+1);
        for i = 1:N
            h_elem = x(i+1) - x(i);
            A_full(i:i+1, i:i+1) = A_full(i:i+1, i:i+1) + [1, -1; -1, 1]/h_elem;
        end
        A = A_full(2:end-1, 2:end-1);

        % Assemble load vector b by mass-lumping
        b_full = zeros(N+1,1);
        for i = 1:N
            h_elem = x(i+1) - x(i);
            b_full(i:i+1) = b_full(i:i+1) + (h_elem/2)*[source(x(i)); source(x(i+1))];
        end
        b = b_full(2:end-1);

        % Solve the FEM system
        u_int = A \ b;
        u_h   = [0; u_int; 0];        % enforce BCs

        % Exact solution and error norms
        u_ex      = sin(pi*x);
        errInf(j) = max(abs(u_ex - u_h));
        w         = u_ex(2:end-1) - u_int;
        errE(j)   = sqrt(w' * A * w);
        errL2(j)  = computeL2(x, u_h);
    end

    % Plot errors
    figure; hold on; box on;
    loglog(H, errL2,  '-o','LineWidth',1.2);
    loglog(H, errE,   '-s','LineWidth',1.2);
    loglog(H, errInf,'-^','LineWidth',1.2);
    legend('L^2 norm','Energy norm','Max norm','Location','SouthEast');
    xlabel('h'); ylabel('Error');
    title('Convergence of FEM for u(x)=sin(\pi x)');
    grid on;

    % Display rates
    displayRates(H, errL2, errE, errInf);
end
```
## Helper Functions
```matlab
function y = source(x)
    % Source term f(x) = π² sin(π x)
    y = pi^2 * sin(pi*x);
end

```


# Usage 
```matlab
% In MATLAB Command Window
FEMPiecewiseLinearErrConv;
```

## License
This project is licensed under the MIT License - see the LICENSE file for details.
```





