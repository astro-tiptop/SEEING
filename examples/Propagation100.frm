Rayleigh-Sommerfeld
Eq(E_1(x_1, y_1, z_1), Integral(-I*E_0*z_1*(I*lambda/(2*pi*sqrt(z_1**2 + (-x_0 + x_1)**2 + (-y_0 + y_1)**2)) + 1)*exp(2*I*pi*sqrt(z_1**2 + (-x_0 + x_1)**2 + (-y_0 + y_1)**2)/lambda)/(lambda*(z_1**2 + (-x_0 + x_1)**2 + (-y_0 + y_1)**2)), (x_0, -a, a), (y_0, -a, a)))
Rayleigh-Sommerfeld Arg
-I*E_0*z_1*(I*lambda/(2*pi*sqrt(z_1**2 + (-x_0 + x_1)**2 + (-y_0 + y_1)**2)) + 1)*exp(2*I*pi*sqrt(z_1**2 + (-x_0 + x_1)**2 + (-y_0 + y_1)**2)/lambda)/(lambda*(z_1**2 + (-x_0 + x_1)**2 + (-y_0 + y_1)**2))
Approximate Rayleigh-Sommerfeld
Eq(E_1(x_1, y_1, z_1), Integral(-I*E_0*z_1*exp(2*I*pi*sqrt(z_1**2 + (-x_0 + x_1)**2 + (-y_0 + y_1)**2)/lambda)/(lambda*(z_1**2 + (-x_0 + x_1)**2 + (-y_0 + y_1)**2)), (x_0, -a, a), (y_0, -a, a)))
Approximate Rayleigh-Sommerfeld Arg
-I*E_0*z_1*exp(2*I*pi*sqrt(z_1**2 + (-x_0 + x_1)**2 + (-y_0 + y_1)**2)/lambda)/(lambda*(z_1**2 + (-x_0 + x_1)**2 + (-y_0 + y_1)**2))
Near Fresnel
Eq(E_1(x_1, y_1, z_1), Integral(-I*E_0*exp(2*I*pi*z_1/lambda)*exp(I*pi*((-x_0 + x_1)**2 + (-y_0 + y_1)**2)/(lambda*z_1))/(lambda*z_1), (x_0, -a, a), (y_0, -a, a)))
Near Fresnel Arg
-I*E_0*exp(2*I*pi*z_1/lambda)*exp(I*pi*((-x_0 + x_1)**2 + (-y_0 + y_1)**2)/(lambda*z_1))/(lambda*z_1)
Far Fresnel
Eq(E_1(x_1, y_1, z_1), Integral(-I*E_0*exp(2*I*pi*z_1/lambda)*exp(I*pi*(x_0**2 + y_0**2)/(lambda*z_1))*exp(2*I*pi*(x_0*x_1 + y_0*y_1)/(lambda*z_1))/(lambda*z_1), (x_0, -a, a), (y_0, -a, a)))
Far Fresnel Arg
-I*E_0*exp(2*I*pi*z_1/lambda)*exp(I*pi*(x_0**2 + y_0**2)/(lambda*z_1))*exp(2*I*pi*(x_0*x_1 + y_0*y_1)/(lambda*z_1))/(lambda*z_1)
Fraunhofer
Eq(E_1(x_1, y_1, z_1), Integral(4*pi**2*E_0*exp(-2*I*pi*(x_0*x_1 + y_0*y_1)/(lambda*z_1))/lambda, (x_0, -a, a), (y_0, -a, a)))
Fraunhofer Arg
4*pi**2*E_0*exp(-2*I*pi*(x_0*x_1 + y_0*y_1)/(lambda*z_1))/lambda
xyCircle
(1 - x_0**2/a**2 - y_0**2/a**2)/(2*Abs(-1 + x_0**2/a**2 + y_0**2/a**2)) + 1/2
xyLens
exp(-I*pi*(x_0**2 + y_0**2)/(2*FN*a*lambda))
