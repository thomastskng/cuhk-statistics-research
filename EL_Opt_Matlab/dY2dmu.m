function result=dY2dmu(t, mu, sigma, delta)
result=-exp(-delta*sigma^2*t^2)*cos(t*(mu-sigma^2/2)*delta)*(t*delta);