#include <nlopt.hpp>
#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* defines */

// comment following line unless compiling for WEB back end
// #define FOR_WEB 1

#define circswap(a,b,c) {double *temp = a; a = b; b = c; c = temp;}
#define swap(a,b)       {double *temp = a; a = b; b = temp;}
#define dmax(a,b) ( (a > b) ? a : b )
#define n1 (n+1)
#define nodeabove(t)  ((int) ( (ceil((t)/dt) > N) ?  N : ceil((t)/dt)))
#define node(t)  ((int) ( (rint((t)/dt) > N) ?  N : rint((t)/dt)))
//#define MAXYLD 100
//#define MAXYRS 100
//#ifdef FOR_WEB
//#define pause() {}
//#endif

// #ifndef FOR_WEB
// #define pause() {if(!file_output){ int ch; \
//                      printf("Press enter to continue ..."); \
// 		     while (EOF != (ch = getchar()) && '\n' != ch);}}
// #endif

/* function declarations */

int allocate(int dim);
double boundary(double st, int n);
void comment(void);
void compute(void);
void dividends(void);
void driver(void);
void forward_rates(void);
int  inconsistent(void);
void init_lattice(void);
void interleave(double yld[], double vol[], int m, double Mat[],
                double Yld[], double Vol[], int NoYields, double Mat2[],
		double Yld2[], double Vol2[],int& k);
double interpolate(double yld[], double t);
void interpolate(double yld[], double dt, int m,
	         double Mat[], double Yld[], int NoYields);
void interpolated_zero_curve(void);
void lattice(void);
void vol_interpolate(double vol[], double dt, int m,
	             double Mat[], double Vol[], int NoVolats);
double node_value(double st, int n, double liveval);
void outputs(void);
void printdata(void);
void print_lattice_yldcve(void);
void print_yldcve(void);
void readdata(void);
void rollback(void);
void set_params(void);
void trace_tree(void);
void yld_volat(double trial, double ratio, double &yield, double &volat);
void zero_curve(double yld[],  double vol[], int m);

/* routines for Embedded Options (EO) */

void EOdriver(void);
double EOnode_value(int n, double liveval);
void EOprintdata(void);
void EOreaddata(void);
void EOset_params(void);

/* global variables */
const int MAXYLD = 100, MAXYRS = 100, nco_max = 20,  npo_max = 20;
double 	*sold, *snew, *s_2, *vold, *vnew, *v_2, *bold, *bnew, *b_2,
  *shigh, *ratios, *yld, *flows, *vol;
//double   *Mat, *Mat2, *Yld, *Yld2, *Vol, *Vol2;
double* Mat = new double[MAXYLD];
double* Mat2 = new double[MAXYLD+2*MAXYRS+1];
double* Yld = new double[MAXYLD];
double* Yld2 = new double[MAXYLD+2*MAXYRS+1];
double* Vol = new double[MAXYLD];
double* Vol2 = new double[MAXYLD+2*MAXYRS+1];
double  redemption_value, straight;
double 	x, Te, t = 0, te, sigma, lvolat_min, yld_error, vol_error, weight, 
  tolerance = 1e-4;
double 	dt, rootdt, p = 0.5;
double 	liveval, delta, delta2, Gamma, theta1, theta2, f, f42,
  fl1, fl2, dsigma, factor0, stoch_durn, stoch_durn2;
double  BondMat, CouponRate, ShortRate;
int 	BondMatDays, BondMatYears, N, N0, n, ne, nTe, opt_type, 
  /*FirstCouponDays,*/ CouponFrequency, YldCurveType, NoYields, NoYields2, 
  minus_one = -1, minT, maxT, nodes;
int 	file_output, put, call, silent = 0, out_data, out_diagnostics,
  out_yldcve, out_lattice_yldcve, out_lattice, out_results, embedded,
  ParBondYldCurve, ParBondVolCurve, fixed_sigma;	/* Boolean */
char    nodechar, *nodechars;
FILE	*infile, *outfile;
// extra variables for embedded options
double 	co_start_0[nco_max], co_end_0[nco_max], co_price[nco_max],
  po_start_0[npo_max], po_end_0[npo_max], po_price[npo_max],
  oas_toler, oas = 0.0, issue_price;
int     co_start[npo_max], co_end[npo_max], po_start[npo_max], po_end[npo_max],
  nco, npo, find_oas, max_iter, iter = 0;

////////////////////////////////////////////////////////////////////
//////////////////   main (driver) routines  ///////////////////////
////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
#ifndef FOR_WEB
  if (argc < 2){
    fprintf(stderr, "Usage BDT InputFileName [OutputFileName]\n");
    exit(1);
  }
  if (NULL == (infile = fopen(argv[1],"r"))){
    fprintf(stderr,"Unable to open %s\n", argv[1]);
    exit(1);
  }
#endif


#ifdef FOR_WEB
  infile = stdin;
  printf("Content-type:text/html\n\n<HEAD>\n<TITLE>%s Valuation Results",
         argv[1]);
  printf("</TITLE></HEAD><BODY>"
         "<P>This software by Prof. J. R. Varma was designed for "
         "educational use and no representation is made that it is "
         "suitable for commercial use. "
         "No warranty is provided that the software is error free.\n"
         "<PRE>\n");
#endif
 
  if (argc < 3){
    outfile = stdout;
    file_output = 0;
  }else{
    file_output = 1;
    if (NULL == (outfile = fopen(argv[2],"w"))){
      fprintf(stderr,"Unable to open %s\n", argv[2]);
      exit(1);
    }
  }
  readdata();
  if (!silent && out_data)
    printdata();
  if(inconsistent())
    exit(1);
  int dim = (NoYields + BondMat*2 > N) ? (int) (NoYields + BondMat*2 + 3): N + 3;
  if (allocate(dim)){
    printf("Out of Memory\n");
    exit(1);
  }
  set_params();
  if (embedded) EOdriver(); else driver();
}

void alloc_char(char* &var, int dim, int& err)
{
  var = (char *) calloc(dim, sizeof(char));
  if(var == NULL)
    err = 1;
  else
    var++;
}

int allocate(int dim)
{
  int err;
#define alloc(var,type) \
	     var = (type *) calloc(dim, sizeof(type)); \
	     if(var == NULL) \
		err = 1; \
	     else     \
		var++;

  err = 0;
  alloc(sold, double);	alloc(snew, double);
  alloc(s_2, double);	alloc(vold, double);
  alloc(vnew, double);	alloc(v_2, double);
  alloc(bold, double);	alloc(bnew, double);
  alloc(b_2, double);	alloc(shigh, double);
  alloc(ratios, double);  alloc(yld, double);
  alloc(flows, double);	alloc(vol, double);
  alloc_char(nodechars,dim,err);
#undef alloc
  return(err);
}

void driver(void)
{
  lattice();
  if (!silent && out_results) outputs();
}


double oas_error(const std::vector<double> &x,
                 std::vector<double> &grad,
                 void *my_func_data) {
  oas = x[0]; 
  driver();
  double error = f-issue_price;
  return sqrt(error * error);
}

void EOdriver()
{
  driver();
  if (!find_oas) return; /* else find OAS */
  double oas_0 = oas;
  nlopt::opt opt(nlopt::LN_NELDERMEAD, 1);
  std::vector<double> x(1);
  double min_oas_error;
  opt.set_min_objective(oas_error, NULL);
  opt.set_ftol_rel(1e-4);
  int silent_0 = silent;
  silent = 1;
  nlopt::result result = opt.optimize(x, min_oas_error);
  silent = silent_0;
  if(result > nlopt::XTOL_REACHED)
    fprintf(outfile, 
	    "OAS computation failed to converge. Results are not valid\\n");
  fprintf(outfile,"Valuation using estimated OAS\n"
	  "OAS = %6.3f%%.  Pricing Error = %6.6G\n"
	  "Detailed results for this value follow.\n\n",
	  oas*100, f-issue_price);
  oas = x[0]; 
  driver();
}

void lattice(void)
{
  init_lattice();
  while(n){
    if (!silent && out_lattice) trace_tree();
    rollback();
  }
  if (!silent && out_lattice) trace_tree();
  compute();
}

void set_params(void)
{
  interpolated_zero_curve();
  forward_rates();
  if (!silent && out_lattice_yldcve) print_lattice_yldcve();
  for(int i = 0; i <= N; i++)
    flows[i] = 0;
  if(CouponRate != 0){
    //    double t = FirstCouponDays/365.0;
    for(double t = BondMat; t > 0; t -= 1.0/CouponFrequency){
      //      printf("t= %f node = %d\n", t, node(t));
      flows[node(t)] += CouponRate/CouponFrequency;
    }
  }
  if (embedded) EOset_params();
}

void EOset_params()
{
  int i;

  for (i = 0; i < nco; i++){
    co_start[i] = node(co_start_0[i]);
    co_end[i] = node(co_end_0[i]);
  }
  for (i = 0; i < npo; i++){
    po_start[i] = node(po_start_0[i]);
    po_end[i] = node(po_end_0[i]);
  }
}

////////////////////////////////////////////////////////////////////
/////////////////   lattice computation routines  //////////////////
////////////////////////////////////////////////////////////////////

void init_lattice(void)
{
  int i;

  n = N;
  for (i = 0; i <= n; i++){
    sold[i] = 0;
    bold[i] = redemption_value;
    vold[i] = (embedded) ? redemption_value :
      (n == nTe)  ? boundary(bold[i], n) : 0;
    nodechars[i] = ' ';
  }
}

void rollback(void)
{
  int i;

  if (n == 4) f42 = vold[2]; /* This value is used for theta estimate */
  n--;
  snew[0] = shigh[n];
  for (i = 1; i <= n; i ++)
    snew[i] = snew[i-1]/ratios[n];
  for (i = 0; i <= n; i++){
    bnew[i] = (p*bold[i] + (1-p)*bold[i+1] + flows[n1]) /
      (1+(snew[i]+oas)*dt);
    double cashflow = (embedded) ? flows[n1] : 0.0;
    liveval = (p*vold[i] + (1-p)*vold[i+1] + cashflow) /
      (1+(snew[i]+oas)*dt);
    vnew[i] = (embedded) ? EOnode_value(n, liveval) :
      node_value(bnew[i], n, liveval);
    nodechars[i] = nodechar;
  }
  circswap(s_2,sold,snew);
  circswap(v_2,vold,vnew);
  circswap(b_2,bold,bnew);
}

double boundary(double bt, int n)
{
  double value;

  if (put)
    value = x - bt;
  else if (call)
    value = bt - x;
  if (n == nTe) value = dmax(0, value);
  return value;
}

double node_value(double bt, int n, double liveval)
{
  double deadval, value;
  int exercise;

  deadval = boundary(bt,n);
  exercise = (n >= ne && n <= nTe && deadval > liveval);
  value = (exercise) ? deadval : liveval;
  nodechar = (exercise) ? 'E' : ' ';
  return(value);
}

double EOnode_value(int n, double liveval)
{
  int call_it, callable, puttable, put_it;

  callable = 0; call_it = -1;
  double callvalue = liveval;
  for (int i = 0; i < nco; i++){
    callable = (n  >= co_start[i]) && (n  <= co_end[i]);
    if (callable && co_price[i] < callvalue){
      call_it = i;
      callvalue = co_price[i];
    }
  }
  puttable = 0; put_it = -1;
  double putvalue = liveval;
  for (int i = 0; i < npo; i++){
    puttable = (n  >= po_start[i]) && (n  <= po_end[i]);
    if (puttable && po_price[i] > putvalue){
      put_it = i;
      putvalue = po_price[i];
    }
  }
  if(put_it >= 0){
    nodechar = 'P' + put_it;
    return putvalue;
  }else if (call_it >= 0){
    nodechar = 'C' + call_it;
    return callvalue;
  }else{
    nodechar = ' ';
    return liveval;
  }
}

void compute(void)
{
  circswap(snew,sold,s_2);
  circswap(bnew,bold,b_2);
  circswap(vnew,vold,v_2); /* reverse the last swap in rollback */
  f = vnew[0] + flows[0];
  if(embedded) straight = bnew[0];
  if(embedded){
    minT = N; maxT = 0;
    for(int i = 0; i < nco; i++){
      if(co_end[i] < minT) minT = co_end[i];
      if(co_end[i] > maxT) maxT = co_end[i];
    }
    for(int i = 0; i < npo; i++){
      if(po_end[i] < minT) minT = po_end[i];
      if(po_end[i] > maxT) maxT = po_end[i];
    }
    nodes = maxT;
  }else
    nodes = nTe;
  if(nodes >= 1) {
    if(!embedded)
      delta = (vold[0] - vold[1])/(bold[0] - bold[1]);
    stoch_durn = -(vold[0] - vold[1])/f/(sold[0] - sold[1]);
  }
  if(nodes >= 2){
    if(!embedded){
      delta2 = (v_2[0] - v_2[2])/(b_2[0] - b_2[2]);
      Gamma = ( (v_2[0] - v_2[1])/(b_2[0] - b_2[1]) -
		(v_2[1] - v_2[2])/(b_2[1] -b_2[2]) ) /
	(0.5 * (b_2[0] - b_2[2]));
    }
    stoch_durn2 = -(v_2[0] - v_2[2])/f/(s_2[0] - s_2[2]);
    theta1 = (v_2[1] + flows[2] - f)/(2*dt);
    if (nodes >= 4) theta2 = (f42 + flows[4] - vnew[0])/(4*dt);
  }
}

////////////////////////////////////////////////////////////////////
////////////////// yield curve related routines  ///////////////////
////////////////////////////////////////////////////////////////////


void interpolated_zero_curve(void)
{
  int i;

  //#define vol_interpolate log_interpolate

  ShortRate *= 0.01;
  oas *= 0.01;
  sigma *= 0.01;
  lvolat_min *= 0.01;
  for(i = 0; i < NoYields; i++){
    Yld[i] *= 0.01;
    if(!fixed_sigma)Vol[i] *= 0.01;
  }
  yld[minus_one] = ShortRate; vol[minus_one] = sigma;
  int np = (int)(2*BondMat+1);
  interpolate(yld, 0.5, np, Mat, Yld, NoYields);
  if(!fixed_sigma)
    vol_interpolate(vol, 0.5, np, Mat, Vol, NoYields);
  zero_curve(yld, vol, np);
  if (!silent && out_yldcve) print_yldcve();
  if(!ParBondYldCurve){
    interpolate(yld, dt, N, Mat, Yld, NoYields);
    vol_interpolate(vol, dt, N, Mat, Vol, NoYields);
  }else{
    interleave(yld, vol, np, Mat, Yld, Vol, NoYields,
	       Mat2, Yld2, Vol2, NoYields2);
    interpolate(yld, dt, N, Mat2, Yld2, NoYields2);
    vol_interpolate(vol, dt, N, Mat2, Vol2, NoYields2);
  }
}
/*
void interpolate(double yld[], double dt, int m,
	         double Mat[], double Yld[], int NoYields)
{
  int i;
  double PrevMat, PrevYld;

  PrevMat = 0;
  n = 0;
  PrevYld = yld[minus_one];
  for(i = 0; i < NoYields; i++){
    for( ; n1*dt <= Mat[i] && n <= m; n++)
      yld[n] = PrevYld + (n1*dt-PrevMat)*(Yld[i]-PrevYld)/(Mat[i]-PrevMat);
    PrevYld = Yld[i];
    PrevMat = Mat[i];
  }
  while(n <= m)
    yld[n++] = PrevYld;
}
void log_interpolate(double vol[], double dt, int m,
	             double Mat[], double Vol[], int NoVolats)
{
  int i;
  double PrevMat, PrevVol, Volat;

  PrevMat = 0;
  n = 0;
  PrevVol = log(vol[minus_one]);
  for(i = 0; i < NoVolats; i++){
    Volat = log(Vol[i]);
    for( ; n1*dt <= Mat[i] && n <= m; n++)
      vol[n] = exp(PrevVol + (n1*dt-PrevMat)*(Volat-PrevVol)/
		   (Mat[i]-PrevMat));
    PrevVol = Volat;
    PrevMat = Mat[i];
  }
  PrevVol = exp(PrevVol);
  while(n <= m)
    vol[n++] = PrevVol;
}
*/

void interpolate(double yld[], double dt, int m,
	         double Mat[], double Yld[], int NoYields)
{
  int i;
  double PrevMat, PrevYld;

  PrevMat = 0;
  n = 0;
  PrevYld = yld[minus_one];
  for(i = 0; i < NoYields; i++){
    double fwd_yld = (Yld[i]*Mat[i] - PrevYld*PrevMat)/(Mat[i]-PrevMat);
    for( ; n1*dt <= Mat[i] && n <= m; n++)
      yld[n] = (PrevYld*PrevMat + fwd_yld*(n1*dt-PrevMat))/(n1*dt);
	//PrevYld + (n1*dt-PrevMat)*(Yld[i]-PrevYld)/(Mat[i]-PrevMat);
    PrevYld = Yld[i];
    PrevMat = Mat[i];
  }
  while(n <= m)
    yld[n++] = PrevYld;
}
void vol_interpolate(double vol[], double dt, int m,
	             double Mat[], double Vol[], int NoVolats)
{
  int i;
  double PrevMat_Sq, PrevVar;

  PrevMat_Sq = 0;
  n = 0;
  PrevVar = vol[minus_one]*vol[minus_one];
  for(i = 0; i < NoVolats; i++){
    double Var = Vol[i]*Vol[i];
    double Mat_Sq = Mat[i]*Mat[i];
    double fwd_var = (Var*Mat_Sq - PrevVar*PrevMat_Sq)/(Mat_Sq - PrevMat_Sq);
    for( ; n1*dt <= Mat[i] && n <= m; n++)
      vol[n] = sqrt(PrevVar*PrevMat_Sq+fwd_var*(n1*dt*n1*dt - PrevMat_Sq))/(n1*dt);
    PrevVar = Var;
    PrevMat_Sq = Mat_Sq;
  }
  double PrevVol = sqrt(PrevVar);
  while(n <= m)
    vol[n++] = PrevVol;
}

void zero_curve(double yld[],  double vol[], int m)
{
  double z, zsum, zhisum, zlosum;
  int i;

  if(!ParBondYldCurve && !ParBondVolCurve)
    return; //there is nothing to do!
  /*
    ri, vi are yields & volatilitilies for zeroes; Ri, Vi for par bonds
    z = zero discount function for n
    Par => Zero
    zsum = sum z.i to n-1
    Rn zsum + (1+Rn) z = 1   (price of par bond = 1)
    => z = (1 - Rn zsum)/(1 + Rn)
    vi = 0.5*log(rhi.i/rlo.i)
    Vi = 0.5*log(Rhi.i/Rlo.i)
    Rhi/Rlo = exp(Vi) and Rhi+Rlo = 2Ri determine Rhi and Rlo as
    Rlo = R/(1+g/2)    and      Rhi = Rlo * (1+g)
    where g = exp(2*Vi) - 1
    We use a weekly lattice for this purpose. Hence the root52 term
    Zero => Par
    zsum = sum z.i to n (not n-1)
    Rn zsum + z = 1   (price of par bond = 1)
    => Rn = (1 - z)/zsum
    vi = 0.5*log(rhi.i/rlo.i)
    Vi = 0.5*log(Rhi.i/Rlo.i)
    rhi/rlo = exp(vi) and rhi+rlo = 2ri determine rhi and rlo as
    rlo = r/(1+g/2)    and      rhi = rlo * (1+g)
    where g = exp(2*vi) - 1
    We use a weekly lattice for this purpose. Hence the root52 term
    The factor of 0.5 is to handle semi-annual yields
  */
  zsum = zhisum = zlosum = 0.0;
  double root52 = sqrt(52.0);
  for(i = 0; i <= m; i++){
    double paryld;
    if(ParBondYldCurve){
      z = (1.0 - 0.5*zsum*yld[i]) / (1.0 + 0.5*yld[i]);
      paryld = yld[i];
      yld[i] = 2.0*(pow(z, -1.0/(i+1)) - 1.0);
      zsum += z;
    }else{
      z = pow(1.0 + 0.5*yld[i], -(i+1));
      zsum += z;
      paryld = 2*(1 - z)/zsum;
    }
    if(ParBondVolCurve){
      double g = exp(2*vol[i]/root52) - 1.0;
      double parlo = paryld/(1.0+0.5*g);
      double parhi = parlo*(1.0+g);
      double zlo = (1.0 - 0.5*zlosum*parlo) / (1.0 + 0.5*parlo);
      double zhi = (1.0 - 0.5*zhisum*parhi) / (1.0 + 0.5*parhi);
      zhisum += zhi;
      zlosum += zlo;
      double zeroHI = 2.0*(pow(zhi, -1.0/(i+1)) - 1.0);
      double zeroLO = 2.0*(pow(zlo, -1.0/(i+1)) - 1.0);
      vol[i] = 0.5*log(zeroHI/zeroLO)*root52;
    }
  }
}

void interleave(double yld[], double vol[], int m, double Mat[],
		double Yld[], double Vol[], int NoYields, double Mat2[],
		double Yld2[], double Vol2[],int& k)
{
  int i;

  i = 0; k = 0;
  while(Mat[i] < 0.5 && i < NoYields){
    Mat2[k] = Mat[i];
    if(!fixed_sigma) Vol2[k] = Vol[i];
    Yld2[k++] = Yld[i++];
  }
  i = 0;
  while(i <= m){
    Mat2[k] = 0.5*(i+1);
    if(!fixed_sigma) Vol2[k] = vol[i];
    Yld2[k++] = yld[i++];
  }
}

static char *buf;
static const int buflen = 1000;
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y)) 

double wtd_error(const std::vector<double> &x,
                   std::vector<double> &grad,
                   void *my_func_data) {
    double yield, volat;
    yld_volat(x[0], x[1], yield, volat);
    double LocVol = 0.5*log(x[1])/sqrt(dt);
    double wt = sqrt(weight);
    double u1 = wt*(yield - yld[n]);
    double u2 = (fixed_sigma)? (LocVol - sigma)/wt : (volat - vol[n])/wt;
    return sqrt(u1*u1 + u2*u2);
}

void forward_rates(void)
{
  double wt, yield, volat, u1, u2, rmax, rmin, tmin, tmax, LocVol, 
    max_yld_err, max_vol_err;

  wt = sqrt(weight);
  buf = new char[buflen];
  bold[minus_one] = bnew[minus_one] = vold[minus_one] = vnew[minus_one] = 0;
  shigh[0] = (pow(1 + 0.5*yld[0], 2*dt) - 1)/dt; // compounding 0.5 => dt
  if (fixed_sigma)
    vol[0] = sigma;
  rmax = exp(2*2*vol[0]*rootdt); // max local volat = 2 volat of short rate
  rmin = exp(2*lvolat_min*rootdt); // user specified min local volatility
  tmin = 1e-5; // spot rate of 0.0 => volatility undefined for n=1
  ratios[0] = exp(2*vol[0]*rootdt);
  factor0 = 1/(1+shigh[0]*dt);
  yld_error = 0.0;
  vol_error = 0.0;
  max_yld_err = 0.0;
  max_vol_err = 0.0;
  if(out_diagnostics){
    fprintf(outfile,"Yield and Volatility Curve Diagnostics\n");
#ifdef pause
    pause();
#endif
  }
  nlopt::opt opt(nlopt::LN_NELDERMEAD, 2);
  std::vector<double> lb(2);
  std::vector<double> ub(2);
  std::vector<double> x(2);
  double min_wtd_error;
  nlopt::result result;
  opt.set_min_objective(wtd_error, NULL);
  opt.set_xtol_abs(tolerance);
  for(n = 1; n < N; n++){
    tmax = 20*shigh[n-1]*rmax;
    x[1]  = ratios[n-1];
    x[0]  = shigh[n-1]*sqrt(x[1]);
    lb[0] = tmin;
    lb[1] = rmin;
    ub[0] = tmax;
    ub[1] = rmax;
    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);
    result = opt.optimize(x, min_wtd_error);
    yld_volat(x[0], x[1], yield, volat);
    LocVol = 0.5*log(x[1])/sqrt(dt);
    shigh[n] = x[0];
    ratios[n] = x[1];
    double yld_err = fabs(yld[n]-yield);
    yld_error += yld_err;
    max_yld_err = MAX(max_yld_err, yld_err);
    double vol_err = (fixed_sigma) ? fabs(LocVol - sigma): 
      fabs(vol[n]-volat);
    vol_error += vol_err;
    max_vol_err = MAX(max_vol_err, vol_err);
    if(out_diagnostics && result > nlopt::XTOL_REACHED)
      fprintf(outfile, "nlopt failed\n"
              "t = %10.6f  YieldError  = %10.3f%%  VolatError = %10.3f%%\n",
              n*dt, 100*yld_err, 100*vol_err);
    if(out_diagnostics && min_wtd_error > tolerance){
      fprintf(outfile,
	      "t = %10.6f  YieldError  = %10.3f%%  VolatError = %10.3f%%\n"
	      "%14s  YldEnvelope = %10.3f%%  LocalVolat = %10.3f%%\n"
	      "%14s  Bounds: [%.3f%%, %.3f%%] and [%.3f%%, %.3f%%]\n",
              n*dt, 100*yld_err, 100*vol_err,
              "", 100*x[0], 100*LocVol,
              "", 100*tmin, 100*tmax, 100*lvolat_min, 200*vol[0]);
#ifdef pause
      pause();
#endif
    }
    yld[n] = yield;
    vol[n] = volat;
    swap(vold, vnew);
    swap(bold, bnew);
  } // for n
  yld_error /= N;
  vol_error /= N;
  if(out_diagnostics || 100*max_yld_err > .01 || 100*max_vol_err > 0.01)
    fprintf(outfile, "%s"
	    "Mean Absolute Error: Yield = %.3f%%. Volatility = %.3f%%\n"
	    "Max  Absolute Error: Yield = %.3f%%. Volatility = %.3f%%\n\n",
	    (100*max_yld_err > .01 || 100*max_vol_err > 0.01) ?
	    "Yield/volatility curve could not be perfectly matched\n" : 
	    "Yield/volatility curve matched satisfactorily\n",
	    100*yld_error, 100*vol_error, 
	    100*max_yld_err, 100*max_vol_err);
}


void yld_volat(double trial, double ratio, double &yield, double &volat)
{
  int i;
  double yldv, yldb, totalv, totalb;

  snew[0] = trial;
  for(i = 1; i <= n; i++)
    snew[i] = snew[i-1] / ratio;
  snew[n+1] = snew[n+2] = totalv = totalb = 0;
  if(n == 1){
    vnew[2] = bnew[0] = 0;
    vnew[0] = vnew[1] = 0.5/(1+snew[0]*dt);
    bnew[1] = bnew[2] = 0.5/(1+snew[1]*dt);
  }
  for(i = 0; i <= n1; i++){
    if(n > 1){
      vnew[i] = 0.5*(vold[i-1]/(1+snew[i-1]*dt) + vold[i]/(1+snew[i]*dt));
      bnew[i] = 0.5*(bold[i-1]/(1+snew[i-1]*dt) + bold[i]/(1+snew[i]*dt));
    }
    totalv += vnew[i];
    totalb += bnew[i];
  }
  yldv = 2*(pow(totalv, -1/(2*(n1-1)*dt)) - 1);
  yldb = 2*(pow(totalb, -1/(2*(n1-1)*dt)) - 1);
  yield = 2*(pow(0.5*(totalv + totalb)*factor0, -1/(2*n1*dt)) - 1);
  if(yldb == 0.0)
    printf("Volatility Undefined\n");
  volat = 0.5 * log(yldv/yldb)/ rootdt;
}


////////////////////////////////////////////////////////////////////
//////////////////   input (read data) routines  ///////////////////
////////////////////////////////////////////////////////////////////

void readdata(void)
{
  int i;

  comment(); fscanf(infile,"%d",&N0);
  N = (N0 < 2) ? 2 : N0;
  comment();
  fscanf(infile,"%d",&out_data);
  fscanf(infile,"%d",&out_yldcve);
  fscanf(infile,"%d",&out_lattice_yldcve);
  fscanf(infile,"%d",&out_diagnostics);
  fscanf(infile,"%d",&out_lattice);
  fscanf(infile,"%d",&out_results);
  // Bond Details
  comment(); fscanf(infile,"%d",&BondMatYears);
  fscanf(infile,"%d",&BondMatDays);
  BondMat = BondMatYears + BondMatDays/365.0;
  dt = BondMat/N;
  rootdt = sqrt(dt);
  comment(); fscanf(infile,"%lf",&CouponRate);
  fscanf(infile,"%d",&CouponFrequency);
  //  fscanf(infile,"%d",&FirstCouponDays);
  comment(); fscanf(infile,"%lf",&redemption_value);
  // Yield Curve
  comment(); fscanf(infile,"%lf",&ShortRate);
  fscanf(infile,"%lf",&sigma);
  comment(); fscanf(infile,"%d",&ParBondYldCurve);
  fscanf(infile,"%d",&fixed_sigma);
  fscanf(infile,"%d",&NoYields);
  ParBondVolCurve = (fixed_sigma == 2);
  fixed_sigma = (fixed_sigma == 0);
  comment();
  for(i = 0; i < NoYields; i++){
    fscanf(infile,"%lf %lf", &Mat[i], &Yld[i]);
    if(!fixed_sigma) fscanf(infile,"%lf", &Vol[i]);
  }
  comment(); fscanf(infile,"%lf",&lvolat_min);
  comment(); fscanf(infile,"%lf",&weight);
  comment(); fscanf(infile,"%lf",&oas);
  // Option Details
  comment(); fscanf(infile,"%d",&opt_type); embedded = (opt_type != 0);
  if (embedded)
    EOreaddata();
  else{
    comment(); fscanf(infile,"%lf",&x);
    comment(); fscanf(infile,"%lf",&Te);   nTe = node(Te);
    comment(); fscanf(infile,"%lf",&te);  ne = node(te);
    comment(); fscanf(infile,"%d",&i); put = (i==0); call = (i==1);
  }
}

void EOreaddata()
{
  int i;

  comment(); fscanf(infile,"%d",&find_oas);
  comment();
  if (find_oas){
    fscanf(infile,"%lf", &issue_price);
  }
  comment(); fscanf(infile,"%d",&nco);
  if (nco > nco_max){
    fprintf(outfile,"Number of Options to Issuer (%d) Exceeds %d\n",
	    nco, nco_max);
    exit(1);
  }
  comment();
  for (i = 0; i < nco; i++)
    fscanf(infile,"%lf %lf %lf",
	   &co_start_0[i], &co_end_0[i], &co_price[i]);
  comment(); fscanf(infile,"%d",&npo);
  if (npo > npo_max){
    fprintf(outfile,"Number of Put Options to Investor (%d) Exceeds %d\n",
	    npo, npo_max);
    exit(1);
  }
  comment();
  for (i = 0; i < npo; i++)
    fscanf(infile,"%lf %lf %lf",
	   &po_start_0[i], &po_end_0[i], &po_price[i]);
}

void comment(void)
{
  int c;

#ifndef newcomment
  while (EOF != (c = getc(infile)) && c != '#');
  while (EOF != (c = getc(infile)) && c != '#');
  if (c == EOF) {
    fprintf(outfile,"End of File During Read\n");
    exit(1);
  }
  //new comment changes above code to:
  //(a) to make the comment optional
  //(b) allow multiple comments one after the other
#endif

#ifdef newcomment
  while (1){
    fscanf(infile," "); c = getc(infile);
    switch(c){
    case EOF:   printf("End of File During Read\n"); exit(1);
    case '#':   while (EOF != (c = getc(infile)) && c != '#'); break;
    default :   ungetc(c, infile); return;
    }
    if (c == EOF){
      printf("End of File During Read\n");
      exit(1);
    }
  }
#endif
}

////////////////////////////////////////////////////////////////////
///////////////////  output (print) routines  //////////////////////
////////////////////////////////////////////////////////////////////

void printdata(void)
{
  int i;

#ifdef pause
  pause();
#endif
  fprintf(outfile,"%20s I n p u t    D a t a\n","");
  // Bond Details
  fprintf(outfile,"\n%20s A. Bond Particulars\n","");
  fprintf(outfile,"%-30s %d years %d days\n",
	  "Bond Maturity", BondMatYears, BondMatDays);
  fprintf(outfile,"%-30s %10.4f%%\n","Coupon Rate ", CouponRate);
  fprintf(outfile,"%-30s %10d\n","Coupon Frequency ", CouponFrequency);
  //  fprintf(outfile,"%-30s %10d\n","Days to 1st Coupon  ", FirstCouponDays);
  fprintf(outfile,"%-30s %10.2f\n","Redemption Value", redemption_value);
  // Yield Curve
  fprintf(outfile,"\n%20s B. Yield Curve\n","");
  fprintf(outfile,"%-30s %10.4f%%\n","Short Rate", ShortRate);
  fprintf(outfile,"%-30s %10.4f%%\n","Short Rate Volatility",sigma);
  fprintf(outfile,"%-30s %10s %10s %10s\n %40s %10s %10s\n", 
	  fixed_sigma ? "Yield Curve" : "Yield/Volatility Curve",
	  "Maturity", "Yield", fixed_sigma ? "" : "Volatility",
	  "", (ParBondYldCurve) ? "(Par Bond)" : "(Zero)",
	  (ParBondVolCurve) ? "(Par Bond)" : fixed_sigma ? "" : "(Zero)");
  for(i = 0; i < NoYields; i++)
    if(fixed_sigma)
      fprintf(outfile,"%30s %10.2f %9.2f%%\n","", Mat[i], Yld[i]);
    else
      fprintf(outfile,"%30s %10.2f %9.2f%% %9.2f%%\n",
	      "", Mat[i], Yld[i], Vol[i]);
  fprintf(outfile,"%-30s %10.4f%%\n","Minimum Local Volatility", lvolat_min);
  fprintf(outfile,"%-30s %10.4f%\n","Weightage (Yield v. Volatility)",
	  weight);
  fprintf(outfile,"%-30s %10.4f%%\n","Spread above Yield Curve", oas);
  // Option Details
  if(!embedded){
    fprintf(outfile,"\n%20sC. Option Particulars\n","");
    fprintf(outfile,"%-30s %10s\n", "Option Type", (put) ? "Put" : "Call");
    fprintf(outfile,"%-30s %10.4f\n","Exercise Price",x);
    fprintf(outfile,"%-30s %10.4f\n","Expires at Time",Te);
    fprintf(outfile,"%-30s %10.4f\n","Exercise at or after Time",te);
  }
  // Computation and Outputs
  fprintf(outfile,"\n%20s D. Computation and Outputs\n","");
  fprintf(outfile,"%-30s %10d\n","No of Lattice Points",N0);
  if (N != N0)
    fprintf(outfile,"%-30s %10d\n","No of Lattice Points Changed to ",N);
  fprintf(outfile,"Outputs:\n");
  if(out_data)
    fprintf(outfile,"%5s%s\n","","Data. ");
  if(out_yldcve)
    fprintf(outfile,"%5s%s\n","","Yield Curve. ");
  switch (out_lattice_yldcve){
  case 0:  break;
  case 1:  fprintf(outfile,"%5s%s\n","","LatticePointWise Yield Curve. "); break;
  default: fprintf(outfile,"%5sLatticePointWise Yield Curve (Every %d points)\n",
		   "", out_lattice_yldcve);
  }
  if(out_diagnostics)
    fprintf(outfile,"%5s%s\n","","Detailed Yield/Volatility Diagnostics ");
  if(out_lattice)
    fprintf(outfile,"%5s%s\n","","Lattice. ");
  if(out_results)
    fprintf(outfile,"%5s%s\n","","Results. ");
  if (embedded) EOprintdata();
  fprintf(outfile,"\n");
}

int inconsistent(void)
{
  int err;
#define error(str) {fprintf(outfile,"Error: %s\n", str); err = 1;}

  err = 0;
  if (BondMat < Te)
    error("Bond Matures before Option Expires");
  if (BondMat < te)
    error("Bond Matures before Option Begins");
  //  if (BondMat < FirstCouponDays/365)
  //  error("Bond Matures before First Coupon");
  if (CouponRate != 0 && CouponFrequency <= 0)
    error("Coupon Frequency Must be Positive for Coupon Bond");
  if (BondMat > MAXYRS){
    fprintf(outfile,"Par Bond to Zero Conversion Can Handle Only %d Years\n",
	    MAXYRS);
    err = 1;
  }
#undef error
  return(err);
}

void EOprintdata()
{
  int i;

  fprintf(outfile,"\n%20s E. Embedded Option Particulars\n","");
  fprintf(outfile,"%-30s %10d\n","No. of call options to issuer",nco);
  for (i = 0; i < nco; i++){
    if (i == 0)
      fprintf(outfile,"%10s %10s %10s\n", "StartTime", "EndTime", "Price");
    fprintf(outfile,"%10.3f %10.3f %10.2f\n", co_start_0[i],
	    co_end_0[i], co_price[i]);
  }
  fprintf(outfile,"%-30s %10d\n","No. of put options to investor",npo);
  for (i = 0; i < npo; i++){
    if (i == 0)
      fprintf(outfile,"%10s %10s %10s\n", "StartTime", "EndTime", "Price");
    fprintf(outfile,"%10.3f %10.3f %10.2f\n", po_start_0[i],
	    po_end_0[i], po_price[i]);
  }
  fprintf(outfile,"%-30s %s\n","Option Adjusted Spread",
          (find_oas ==0) ? "Not Required" : "Required");
  if(find_oas){
    fprintf(outfile,"%-30s %10.2f\n","Issue Price", issue_price);
  }
}

void print_yldcve(void)
{
#ifdef pause
  pause();
#endif
  fprintf(outfile,"%20s Y i  e l d   C u r v e\n\n","");
  if(!fixed_sigma)
    fprintf(outfile,"%10s %10s %10s %10s %10s\n",
	    "Maturity", "ZeroYield", "ZeroVolat", "ParYield", "ParVolat");
  else
    fprintf(outfile,"%10s %10s %10s\n",
	    "Maturity", "ZeroYield", "ParYield");
  double zsum = 0; double zhisum = 0;  double zlosum = 0;
  double root52 = sqrt(52.0);
  for(int t = 0; t < 2*BondMat; t ++){
    /*
      ri, vi are yields & volatilitilies for zeroes; Ri, Vi for par bonds
      z = zero discount function for n
      zsum = sum z.i to n (not n-1)
      Rn zsum + z = 1   (price of par bond = 1)
      => Rn = (1 - z)/zsum
      vi = 0.5*log(rhi.i/rlo.i)
      Vi = 0.5*log(Rhi.i/Rlo.i)
      rhi/rlo = exp(vi) and rhi+rlo = 2ri determine rhi and rlo as
      rlo = r/(1+g/2)    and      rhi = rlo * (1+g)
      where g = exp(2*vi) - 1
      We use a weekly lattice for this purpose. Hence the root52 term
    */
    double t1 = t+1;
    double z = 1.0 /pow(1+0.5*yld[t],t1); // 0.5*yld = semi-annual coup
    zsum += z;
    double paryld = 2*(1 - z)/zsum;
    if(!fixed_sigma){
      double g = exp(2*vol[t]/root52) - 1.0;
      double zlo = 1.0 /pow(1+0.5*yld[t]/(1.0+0.5*g),t1);
      double zhi = 1.0 /pow(1+0.5*yld[t]*(1.0+g)/(1.0+0.5*g),t1);
      zhisum += zhi;
      double paryldHI = 2*(1 - zhi)/zhisum;
      zlosum += zlo;
      double paryldLO = 2*(1 - zlo)/zlosum;
      double parvol = 0.5*log(paryldHI/paryldLO)*root52;
      fprintf(outfile,"%10.2f %9.2f%% %9.2f%% %9.2f%% %9.2f%%\n",
	      0.5*t1, 100*yld[t], 100*vol[t], 100*paryld, 100*parvol);
    }else
      fprintf(outfile,"%10.2f %9.2f%% %9.2f%%\n",
	      t1*dt, 100*yld[t], 100*paryld);
  }// for t
  fprintf(outfile,"\n");
}

void print_lattice_yldcve(void)
{
#ifdef pause
  pause();
#endif
  fprintf(outfile,"%20s Lattice-point-wise Yield Curve\n\n","");
  fprintf(outfile,"%10s %10s %10s %10s %10s\n",
	  "Maturity", "ZeroYield", "FwdRate", "ZeroVolat", "LocalVolat");
  double zlag = 1;
  double rootdt = sqrt(dt);
  for(int t = 0; t < N; t ++){
    double localvolat = 0.5*log(ratios[t])/rootdt;
    double z = pow(1 + yld[t]*dt, -(t+1));
    double fwdrate = (zlag/z - 1)/dt;
    zlag = z;
    if(t % out_lattice_yldcve == 0 || t == N-1)
      fprintf(outfile,"%10.2f %9.2f%% %9.2f%% %9.2f%% %9.2f%%\n",
	      (t+1)*dt, 100*yld[t], 100*fwdrate, 100*vol[t],
	      100*localvolat);
  }// for t
  fprintf(outfile,"\n");
}

double interpolate(double yld[], double t)
{
  int above = nodeabove(t); 
  int below = above - 1;
  if(above > N)
    return(yld[below]);
  else if(below < 0)
    return(yld[above]);
  else
    return(yld[below] + (yld[above] - yld[below])*(t/dt-below));
}

void trace_tree(void)
{
  int i;

  if(n == N){
#ifdef pause
    pause();
#endif
    fprintf(outfile,"%25s L a t t i c e\n\n","");
  }
  if(!embedded)
    fprintf(outfile,"Step = %d (Time = %.3f). Legend: [Rate, Bond, Option]\n",
	    n,n*dt);
  else
    fprintf(outfile,"Step = %d (Time = %.3f). Legend: [Rate, Straight, Bond]\n",
	    n,n*dt);
  fprintf(outfile,"Cash Flow = %.3f\n", flows[n]);
  for (i = 0; i <= n; i++){
    fprintf(outfile,"[%.2f%%, %.4G, %.4G]%c ",
	    sold[i]*100, bold[i], vold[i], nodechars[i]);
    if ( (i % 4 == 3) || (i == n) ) fprintf(outfile,"\n");
  }
  fprintf(outfile,"\n");
}

void outputs(void)
{
#ifdef pause
  pause();
#endif
  fprintf(outfile,"%10s V a l u a t i o n \n\n","");
  fprintf(outfile,"%-30s = %9.4f%%\n","Option Adjusted Spread", 100 * oas);
  if (embedded){
    fprintf(outfile,"%-30s = %10.4f\n","Bond Value with Calls/Puts", f);
    fprintf(outfile,"%-30s = %10.4f\n","Value of Straight Bond", 
	    straight+flows[0]);
  }else
    fprintf(outfile,"%-30s = %10.4f\n","Option Value", f);
  if (nodes >= 1 && !embedded)
    fprintf(outfile,"%-30s = %10.4f(time %.6G to time %.6G)\n",
	    "Delta", delta, t, t + dt);
  if (nodes >= 2 && !embedded)
    fprintf(outfile,"%-30s = %10.4f(time %.6G to time %.6G)\n",
	    "", delta2, t, t + 2*dt);
  if (nodes >= 1)
    fprintf(outfile,"%-30s = %10.4f(using time %.6G to time %.6G)\n",
	    "Stochastic Duration", stoch_durn, t, t + dt);
  if (nodes >= 2)
    fprintf(outfile,"%-30s = %10.4f(using time %.6G to time %.6G)\n",
	    "", stoch_durn2, t, t + 2*dt);
  if (nodes >= 2) {
    if(!embedded)
      fprintf(outfile,"%-30s = %10.4f(time %.6G to time %.6G)\n",
	      "Gamma", Gamma, t, t + 2*dt);
    fprintf(outfile,"%-30s = %10.4f(%.2f%%)\n"
	    "%-30s   (using time %.6G to time %.6G)\n",
	    "Theta", theta1, 100*theta1/f, "", t, t + 2*dt);
  }
  if (nodes >=4)
    fprintf(outfile,"%-30s = %10.4f(%.2f%%)\n"
	    "%-30s   (using time %.6G to time %.6G)\n",
	    "Theta", theta2, 100*theta2/f, "", t, t + 4*dt);
  if (nodes < 2)
    fprintf(outfile,"%s.  %s%sGamma and Theta not estimated\n",
	    (embedded && nco+npo == 0) ? "No embedded options" : 
	    "Lattice too coarse",
	    (nodes < 1) ? "Stochastic Duration, " : "",
	    (nodes < 1 && !embedded) ? "Delta, " : "");
  if (embedded && (minT < 4 && nodes >= 4 || minT < 2 && nodes >= 2) ){
    fprintf(outfile,"Since lattice was coarse, the effect of some of the "
	    "embedded options \nin the bond are not reflected in the "
	    "above estimate of ");
    if(minT < 2 && nodes >=4)
      fprintf(outfile,"theta (both estimates)\n");
    if(minT < 2 && nodes >=2  && nodes < 4)
      fprintf(outfile,"theta\n");
    if(minT >=2 && minT < 4 && nodes >=4)
      fprintf(outfile,"theta (2nd estimate)\n");
  }
}


