pCorr<-function(x,y,prob){
mu_x = sum(x*prob);
mu_y = sum(y*prob);
mu_xy = sum(x*y*prob);
cov_xy = mu_xy - mu_x*mu_y;

var_x = sum((x^2)*prob) - mu_x^2;
var_y = sum((y^2)*prob) - mu_y^2;
corr = cov_xy/(sqrt(var_x*var_y));
return(corr);
}