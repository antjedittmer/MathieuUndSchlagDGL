function dA_dpsi = myJacobian(psi, x, gamma, mu_param, d2, d3, ebeta)
    % Derivative of a11 with respect to psi
    da11_dpsi = (gamma / 2) * mu_param * d3 * cos(psi);

    % Derivative of a12 with respect to psi
    da12_dpsi = (gamma / 2) * (-mu_param * (d3 + ebeta * d2) * sin(psi) + mu_param^2 * d2 * cos(2 * psi));

    % Derivative of a21 and a22 with respect to psi (both are constants)
    da21_dpsi = 0;
    da22_dpsi = 0;

    % Construct the Jacobian matrix
    dA_dpsi = -[da11_dpsi, da12_dpsi;
                da21_dpsi, da22_dpsi];
end