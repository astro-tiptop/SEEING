noisePropagationCoefficient
Eq(C, 0.148843545191282*pi**2*D**2/N_sa_tot)
noisePSDTip
Eq(phi^noise_Tip^f(f), C*sigma^2_WCoG/(Delta_F*df*mu_WCoG**2))
noisePSDTilt
Eq(phi^noise_Tilt^f(f), C*sigma^2_WCoG/(Delta_F*df*mu_WCoG**2))
turbPSDTip
Eq(phi^turb_Tip^f(f), Integral(0.3664*f**2*besselj(2, 2*pi*R*sqrt(k_y**2 + f**2/V**2))**2/(pi**2*R**2*V**3*r_0**(5/3)*(k_y**2 + f**2/V**2)**2*((1.0/L_0)**2 + k_y**2 + f**2/V**2)**(11/6)), (k_y, k_y_min, k_y_max)))
turbPSDTilt
Eq(phi^turb_Tilt^f(f), Integral(0.3664*(1.0 - f**2/(V**2*(k_y**2 + f**2/V**2)))*besselj(2, 2*pi*R*sqrt(k_y**2 + f**2/V**2))**2/(pi**2*R**2*V*r_0**(5/3)*(k_y**2 + f**2/V**2)*((1.0/L_0)**2 + k_y**2 + f**2/V**2)**(11/6)), (k_y, k_y_min, k_y_max)))
interactionMatrixNGS
Matrix([[1, 0, 4*sqrt(3)*D*H_DM*x_NGS/(D + 2*H_DM*r_FoV)**2, 2*sqrt(6)*D*H_DM*y_NGS/(D + 2*H_DM*r_FoV)**2, 2*sqrt(6)*D*H_DM*x_NGS/(D + 2*H_DM*r_FoV)**2], [0, 1, 4*sqrt(3)*D*H_DM*y_NGS/(D + 2*H_DM*r_FoV)**2, 2*sqrt(6)*D*H_DM*x_NGS/(D + 2*H_DM*r_FoV)**2, -2*sqrt(6)*D*H_DM*y_NGS/(D + 2*H_DM*r_FoV)**2]])
residualTT
Eq(res, sqrt(epsilon_Tilt + epsilon_Tip))
residualTip
Eq(epsilon_Tip, Integral(phi^res_Tip, (f, f_min, f_max)))
residualTilt
Eq(epsilon_Tilt, Integral(phi^res_Tilt, (f, f_min, f_max)))
residualTipPSD
Eq(phi^res_Tip, phi^noise_Tip*Abs(H^N_Tip)**2 + phi^wind_Tip*Abs(H^R_Tip)**2)
residualTiltPSD
Eq(phi^res_Tilt, phi^noise_Tilt*Abs(H^N_Tilt)**2 + phi^wind_Tilt*Abs(H^R_Tilt)**2)
ztfTipWindMono
Eq(H^R_Tip^z(z), (1 - 1/z)/(g^Tip_0*z**(-d) + 1 - 1/z))
ztfTiltWindMono
Eq(H^R_Tilt^z(z), (1 - 1/z)/(g^Tilt_0*z**(-d) + 1 - 1/z))
ztfTipNoiseMono
Eq(H^N_Tip^z(z), g^Tip_0*z**(-d)/(g^Tip_0*z**(-d) + 1 - 1/z))
ztfTiltNoiseMono
Eq(H^N_Tilt^z(z), g^Tilt_0*z**(-d)/(g^Tilt_0*z**(-d) + 1 - 1/z))
ztfTipWind
Eq(H^R_Tip^z(z), (1 - 1/z)**2/((g^Tip_1 + 1 - 1/z)*(g^Tip_0*z**(-d) + 1 - 1/z)))
ztfTiltWind
Eq(H^R_Tilt^z(z), (1 - 1/z)**2/((g^Tilt_1 + 1 - 1/z)*(g^Tilt_0*z**(-d) + 1 - 1/z)))
ztfTipNoise
Eq(H^N_Tip^z(z), g^Tip_0*g^Tip_1*z**(-d)/((g^Tip_1 + 1 - 1/z)*(g^Tip_0*z**(-d) + 1 - 1/z)))
ztfTiltNoise
Eq(H^N_Tilt^z(z), g^Tilt_0*g^Tilt_1*z**(-d)/((g^Tilt_1 + 1 - 1/z)*(g^Tilt_0*z**(-d) + 1 - 1/z)))
tfTipWind
Eq(H^R_Tip^f(f), (1 - exp(-2*I*pi*f/f_loop))**2/((g^Tip_1 + 1 - exp(-2*I*pi*f/f_loop))*(g^Tip_0*exp(-2*I*pi*d*f/f_loop) + 1 - exp(-2*I*pi*f/f_loop))))
tfTiltWind
Eq(H^R_Tilt^f(f), (1 - exp(-2*I*pi*f/f_loop))**2/((g^Tilt_1 + 1 - exp(-2*I*pi*f/f_loop))*(g^Tilt_0*exp(-2*I*pi*d*f/f_loop) + 1 - exp(-2*I*pi*f/f_loop))))
tfTipNoise
Eq(H^N_Tip^f(f), g^Tip_0*g^Tip_1*exp(-2*I*pi*d*f/f_loop)/((g^Tip_1 + 1 - exp(-2*I*pi*f/f_loop))*(g^Tip_0*exp(-2*I*pi*d*f/f_loop) + 1 - exp(-2*I*pi*f/f_loop))))
tfTiltNoise
Eq(H^N_Tilt^f(f), g^Tilt_0*g^Tilt_1*exp(-2*I*pi*d*f/f_loop)/((g^Tilt_1 + 1 - exp(-2*I*pi*f/f_loop))*(g^Tilt_0*exp(-2*I*pi*d*f/f_loop) + 1 - exp(-2*I*pi*f/f_loop))))
completeIntegralTipLO
Integral(g^Tip_0**2*phi^noise_Tip/(g^Tip_0**2 - g^Tip_0*exp(2*I*pi*f/f_loop)*exp(-2*I*pi*d*f/f_loop) + g^Tip_0*exp(2*I*pi*d*f/f_loop) + g^Tip_0*exp(-2*I*pi*d*f/f_loop) - g^Tip_0*exp(-2*I*pi*f/f_loop)*exp(2*I*pi*d*f/f_loop) - exp(2*I*pi*f/f_loop) + 2 - exp(-2*I*pi*f/f_loop)) + phi^wind_Tip*(-exp(2*I*pi*f/f_loop) + 2 - exp(-2*I*pi*f/f_loop))/(g^Tip_0**2 - g^Tip_0*exp(2*I*pi*f/f_loop)*exp(-2*I*pi*d*f/f_loop) + g^Tip_0*exp(2*I*pi*d*f/f_loop) + g^Tip_0*exp(-2*I*pi*d*f/f_loop) - g^Tip_0*exp(-2*I*pi*f/f_loop)*exp(2*I*pi*d*f/f_loop) - exp(2*I*pi*f/f_loop) + 2 - exp(-2*I*pi*f/f_loop)), (f, f_min, f_max))
completeIntegralTiltLO
Integral(g^Tilt_0**2*phi^noise_Tilt/(g^Tilt_0**2 - g^Tilt_0*exp(2*I*pi*f/f_loop)*exp(-2*I*pi*d*f/f_loop) + g^Tilt_0*exp(2*I*pi*d*f/f_loop) + g^Tilt_0*exp(-2*I*pi*d*f/f_loop) - g^Tilt_0*exp(-2*I*pi*f/f_loop)*exp(2*I*pi*d*f/f_loop) - exp(2*I*pi*f/f_loop) + 2 - exp(-2*I*pi*f/f_loop)) + phi^wind_Tilt*(-exp(2*I*pi*f/f_loop) + 2 - exp(-2*I*pi*f/f_loop))/(g^Tilt_0**2 - g^Tilt_0*exp(2*I*pi*f/f_loop)*exp(-2*I*pi*d*f/f_loop) + g^Tilt_0*exp(2*I*pi*d*f/f_loop) + g^Tilt_0*exp(-2*I*pi*d*f/f_loop) - g^Tilt_0*exp(-2*I*pi*f/f_loop)*exp(2*I*pi*d*f/f_loop) - exp(2*I*pi*f/f_loop) + 2 - exp(-2*I*pi*f/f_loop)), (f, f_min, f_max))
completeIntegralTip
Integral(g^Tip_0**2*g^Tip_1**2*phi^noise_Tip/((g^Tip_1**2 - g^Tip_1*exp(2*I*pi*f/f_loop) + 2*g^Tip_1 - g^Tip_1*exp(-2*I*pi*f/f_loop) - exp(2*I*pi*f/f_loop) + 2 - exp(-2*I*pi*f/f_loop))*(g^Tip_0**2 - g^Tip_0*exp(2*I*pi*f/f_loop)*exp(-2*I*pi*d*f/f_loop) + g^Tip_0*exp(2*I*pi*d*f/f_loop) + g^Tip_0*exp(-2*I*pi*d*f/f_loop) - g^Tip_0*exp(-2*I*pi*f/f_loop)*exp(2*I*pi*d*f/f_loop) - exp(2*I*pi*f/f_loop) + 2 - exp(-2*I*pi*f/f_loop))) + phi^wind_Tip*(1 - exp(-2*I*pi*f/f_loop))**2*(exp(2*I*pi*f/f_loop) - 1)**2/((g^Tip_1**2 - g^Tip_1*exp(2*I*pi*f/f_loop) + 2*g^Tip_1 - g^Tip_1*exp(-2*I*pi*f/f_loop) - exp(2*I*pi*f/f_loop) + 2 - exp(-2*I*pi*f/f_loop))*(g^Tip_0**2 - g^Tip_0*exp(2*I*pi*f/f_loop)*exp(-2*I*pi*d*f/f_loop) + g^Tip_0*exp(2*I*pi*d*f/f_loop) + g^Tip_0*exp(-2*I*pi*d*f/f_loop) - g^Tip_0*exp(-2*I*pi*f/f_loop)*exp(2*I*pi*d*f/f_loop) - exp(2*I*pi*f/f_loop) + 2 - exp(-2*I*pi*f/f_loop))), (f, f_min, f_max))
completeIntegralTilt
Integral(g^Tilt_0**2*g^Tilt_1**2*phi^noise_Tilt/((g^Tilt_1**2 - g^Tilt_1*exp(2*I*pi*f/f_loop) + 2*g^Tilt_1 - g^Tilt_1*exp(-2*I*pi*f/f_loop) - exp(2*I*pi*f/f_loop) + 2 - exp(-2*I*pi*f/f_loop))*(g^Tilt_0**2 - g^Tilt_0*exp(2*I*pi*f/f_loop)*exp(-2*I*pi*d*f/f_loop) + g^Tilt_0*exp(2*I*pi*d*f/f_loop) + g^Tilt_0*exp(-2*I*pi*d*f/f_loop) - g^Tilt_0*exp(-2*I*pi*f/f_loop)*exp(2*I*pi*d*f/f_loop) - exp(2*I*pi*f/f_loop) + 2 - exp(-2*I*pi*f/f_loop))) + phi^wind_Tilt*(1 - exp(-2*I*pi*f/f_loop))**2*(exp(2*I*pi*f/f_loop) - 1)**2/((g^Tilt_1**2 - g^Tilt_1*exp(2*I*pi*f/f_loop) + 2*g^Tilt_1 - g^Tilt_1*exp(-2*I*pi*f/f_loop) - exp(2*I*pi*f/f_loop) + 2 - exp(-2*I*pi*f/f_loop))*(g^Tilt_0**2 - g^Tilt_0*exp(2*I*pi*f/f_loop)*exp(-2*I*pi*d*f/f_loop) + g^Tilt_0*exp(2*I*pi*d*f/f_loop) + g^Tilt_0*exp(-2*I*pi*d*f/f_loop) - g^Tilt_0*exp(-2*I*pi*f/f_loop)*exp(2*I*pi*d*f/f_loop) - exp(2*I*pi*f/f_loop) + 2 - exp(-2*I*pi*f/f_loop))), (f, f_min, f_max))
ZernikeCovarianceD
Eq(dW_phi(rho), 0.0229*(-1)**m_k*2**(-0.5*KroneckerDelta(0, m_j) - 0.5*KroneckerDelta(0, m_k) + 1)*I**(n_j + n_k)*sqrt((n_j + 1)*(n_k + 1))*(I**(3*m_j + 3*m_k)*cos(theta*(m_j + m_k) + pi*((1 - KroneckerDelta(0, m_j))*((-1)**j - 1) + (1 - KroneckerDelta(0, m_k))*((-1)**k - 1))/4)*besselj(m_j + m_k, 2*pi*f*h*rho) + I**(3*Abs(m_j - m_k))*cos(theta*(m_j - m_k) + pi*((1 - KroneckerDelta(0, m_j))*((-1)**j - 1) - (1 - KroneckerDelta(0, m_k))*((-1)**k - 1))/4)*besselj(Abs(m_j - m_k), 2*pi*f*h*rho))*besselj(n_j + 1, 2*pi*R_1*f)*besselj(n_j + 1, 2*pi*R_2*f)/(pi*R_1*R_2*f*r_0**(5/3)*((1/L_0)**2 + f**2)**(11/6)))
ZernikeCovarianceI
Eq(W_phi(rho), Integral(0.0229*(-1)**m_k*2**(-0.5*KroneckerDelta(0, m_j) - 0.5*KroneckerDelta(0, m_k) + 1)*I**(n_j + n_k)*sqrt((n_j + 1)*(n_k + 1))*(I**(3*m_j + 3*m_k)*cos(theta*(m_j + m_k) + pi*((1 - KroneckerDelta(0, m_j))*((-1)**j - 1) + (1 - KroneckerDelta(0, m_k))*((-1)**k - 1))/4)*besselj(m_j + m_k, 2*pi*f*h*rho) + I**(3*Abs(m_j - m_k))*cos(theta*(m_j - m_k) + pi*((1 - KroneckerDelta(0, m_j))*((-1)**j - 1) - (1 - KroneckerDelta(0, m_k))*((-1)**k - 1))/4)*besselj(Abs(m_j - m_k), 2*pi*f*h*rho))*besselj(n_j + 1, 2*pi*R_1*f)*besselj(n_j + 1, 2*pi*R_2*f)/(pi*R_1*R_2*f*r_0**(5/3)*((1/L_0)**2 + f**2)**(11/6)), (f, f_min, f_max)))
TruncatedMeanBasic
Eq(mu_k_thr, Sum((b + f_k)**i*(nu*(erf(sqrt(2)*(b - i + t)/(2*sigma_RON))/2 + 1/2) + sqrt(2)*sigma_RON*exp(-(b - i + t)**2/(2*sigma_RON**2))/(2*sqrt(pi)) + (-b + i)*(erf(sqrt(2)*(-b + i - t)/(2*sigma_RON))/2 + 1/2))*exp(-b - f_k)/factorial(i), (i, 0, i_max)))
TruncatedVarianceBasic
Eq(sigma^2_k_thr, -mu_k_thr**2 + Sum((b + f_k)**i*(nu**2*(erf(sqrt(2)*(b - i + t)/(2*sigma_RON))/2 + 1/2) + sqrt(2)*sigma_RON*(-b + i + t)*exp(-(b - i + t)**2/(2*sigma_RON**2))/(2*sqrt(pi)) + (sigma_RON**2 + (-b + i)**2)*(erf(sqrt(2)*(-b + i - t)/(2*sigma_RON))/2 + 1/2))*exp(-b - f_k)/factorial(i), (i, 0, i_max)))
TruncatedMean
Eq(mu_k_thr, (-b*(1/2 - erf(-sqrt(2)*(-b - t)/(2*sigma_RON))/2) + nu*(erf(sqrt(2)*(b + t)/(2*sigma_RON))/2 + 1/2) + sqrt(2)*sigma_RON*exp(-(b + t)**2/(2*sigma_RON**2))/(2*sqrt(pi)))*exp(-b - f_k) + Sum((b + f_k)**i*exp(-b - f_k)*Integral(z**(i/(F - 1) - 1)*(nu*(erf(sqrt(2)*(b + t - z)/(2*sigma_RON))/2 + 1/2) + sqrt(2)*sigma_RON*exp(-(b + t - z)**2/(2*sigma_RON**2))/(2*sqrt(pi)) + (-b + z)*(erf(sqrt(2)*(-b - t + z)/(2*sigma_RON))/2 + 1/2))*exp(-z/(F - 1))*exp(-i*log(F - 1)/(F - 1))/gamma(i/(F - 1)), (z, 0.0001, z_max))/factorial(i), (i, 1, i_max)))
TruncatedMeanIntegrand
Eq(I_mu_k_thr, z**(i/(F - 1) - 1)*(nu*(erf(sqrt(2)*(b + t - z)/(2*sigma_RON))/2 + 1/2) + sqrt(2)*sigma_RON*exp(-(b + t - z)**2/(2*sigma_RON**2))/(2*sqrt(pi)) + (-b + z)*(erf(sqrt(2)*(-b - t + z)/(2*sigma_RON))/2 + 1/2))*exp(-z/(F - 1))*exp(-i*log(F - 1)/(F - 1))/gamma(i/(F - 1)))
TruncatedVariance
Eq(sigma^2_k_thr, -mu_k_thr**2 + (nu**2*(erf(sqrt(2)*(b + t)/(2*sigma_RON))/2 + 1/2) + sqrt(2)*sigma_RON*(-b + t)*exp(-(b + t)**2/(2*sigma_RON**2))/(2*sqrt(pi)) + (1/2 - erf(-sqrt(2)*(-b - t)/(2*sigma_RON))/2)*(b**2 + sigma_RON**2))*exp(-b - f_k) + Sum((b + f_k)**i*exp(-b - f_k)*Integral(z**(i/(F - 1) - 1)*(nu**2*(erf(sqrt(2)*(b + t - z)/(2*sigma_RON))/2 + 1/2) + sqrt(2)*sigma_RON*(-b + t + z)*exp(-(b + t - z)**2/(2*sigma_RON**2))/(2*sqrt(pi)) + (sigma_RON**2 + (-b + z)**2)*(erf(sqrt(2)*(-b - t + z)/(2*sigma_RON))/2 + 1/2))*exp(-z/(F - 1))*exp(-i*log(F - 1)/(F - 1))/gamma(i/(F - 1)), (z, 0.0001, z_max))/factorial(i), (i, 1, i_max)))
TruncatedVarianceIntegrand
Eq(I_sigma_k_thr, z**(i/(F - 1) - 1)*(nu**2*(erf(sqrt(2)*(b + t - z)/(2*sigma_RON))/2 + 1/2) + sqrt(2)*sigma_RON*(-b + t + z)*exp(-(b + t - z)**2/(2*sigma_RON**2))/(2*sqrt(pi)) + (sigma_RON**2 + (-b + z)**2)*(erf(sqrt(2)*(-b - t + z)/(2*sigma_RON))/2 + 1/2))*exp(-z/(F - 1))*exp(-i*log(F - 1)/(F - 1))/gamma(i/(F - 1)))
truncatedMeanComponents0
(-b*(1/2 - erf(-sqrt(2)*(-b - t)/(2*sigma_RON))/2) + nu*(erf(sqrt(2)*(b + t)/(2*sigma_RON))/2 + 1/2) + sqrt(2)*sigma_RON*exp(-(b + t)**2/(2*sigma_RON**2))/(2*sqrt(pi)))*exp(-b - f_k)
truncatedMeanComponents1
(b + f_k)**i*exp(-b - f_k)/factorial(i)
truncatedMeanComponents2
Integral(z**(i/(F - 1) - 1)*(nu*(erf(sqrt(2)*(b + t - z)/(2*sigma_RON))/2 + 1/2) + sqrt(2)*sigma_RON*exp(-(b + t - z)**2/(2*sigma_RON**2))/(2*sqrt(pi)) + (-b + z)*(erf(sqrt(2)*(-b - t + z)/(2*sigma_RON))/2 + 1/2))*exp(-z/(F - 1))*exp(-i*log(F - 1)/(F - 1))/gamma(i/(F - 1)), (z, 0.0001, z_max))
truncatedVarianceComponents0
(nu**2*(erf(sqrt(2)*(b + t)/(2*sigma_RON))/2 + 1/2) + sqrt(2)*sigma_RON*(-b + t)*exp(-(b + t)**2/(2*sigma_RON**2))/(2*sqrt(pi)) + (1/2 - erf(-sqrt(2)*(-b - t)/(2*sigma_RON))/2)*(b**2 + sigma_RON**2))*exp(-b - f_k)
truncatedVarianceComponents1
(b + f_k)**i*exp(-b - f_k)/factorial(i)
truncatedVarianceComponents2
Integral(z**(i/(F - 1) - 1)*(nu**2*(erf(sqrt(2)*(b + t - z)/(2*sigma_RON))/2 + 1/2) + sqrt(2)*sigma_RON*(-b + t + z)*exp(-(b + t - z)**2/(2*sigma_RON**2))/(2*sqrt(pi)) + (sigma_RON**2 + (-b + z)**2)*(erf(sqrt(2)*(-b - t + z)/(2*sigma_RON))/2 + 1/2))*exp(-z/(F - 1))*exp(-i*log(F - 1)/(F - 1))/gamma(i/(F - 1)), (z, 0.0001, z_max))
