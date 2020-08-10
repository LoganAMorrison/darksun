#ifndef DARKSUN_STIFF_RADAU_HPP
#define DARKSUN_STIFF_RADAU_HPP

#include "darksun/stiff/common.hpp"
#include "darksun/stiff/decomp.hpp"
#include "darksun/stiff/decsol.hpp"
#include "darksun/stiff/radau_error_estimate.hpp"
#include "darksun/stiff/radau_solvers.hpp"
#include "darksun/stiff/vector_matrix.hpp"
#include <cmath>

namespace darksun::stiff {

struct RadauWeight {
  int nn, ns;
  double xsol, hsol, c[8];
};

struct RadauCoeff3 {
  double t311, t312, t313, t321, t322, t323, t331, ti311, ti312, ti313, ti321,
      ti322, ti323, ti331, ti332, ti333;
};

struct RadauCoeff5 {
  double t511, t512, t513, t514, t515, t521, t522, t523, t524, t525, t531, t532,
      t533, t534, t535, t541, t542, t543, t544, t545, t551, ti511, ti512, ti513,
      ti514, ti515, ti521, ti522, ti523, ti524, ti525, ti531, ti532, ti533,
      ti534, ti535, ti541, ti542, ti543, ti544, ti545, ti551, ti552, ti553,
      ti554, ti555;
};

struct RadauCoeff7 {
  double t711, t712, t713, t714, t715, t716, t717, t721, t722, t723, t724, t725,
      t726, t727, t731, t732, t733, t734, t735, t736, t737, t741, t742, t743,
      t744, t745, t746, t747, t751, t752, t753, t754, t755, t756, t757, t761,
      t762, t763, t764, t765, t766, t767, t771, ti711, ti712, ti713, ti714,
      ti715, ti716, ti717, ti721, ti722, ti723, ti724, ti725, ti726, ti727,
      ti731, ti732, ti733, ti734, ti735, ti736, ti737, ti741, ti742, ti743,
      ti744, ti745, ti746, ti747, ti751, ti752, ti753, ti754, ti755, ti756,
      ti757, ti761, ti762, ti763, ti764, ti765, ti766, ti767, ti771, ti772,
      ti773, ti774, ti775, ti776, ti777;
};

/* Table of constant values */

int c__9 = 9;
int c__1 = 1;
int c__3 = 3;
int c__5 = 5;
double c_b90 = 1.;
double c_b100 = .8;
double c_b140 = 9.;
double c_b141 = .33333333333333331;

void coertv(const int nsmax, RadauCoeff3 &rcoeff3, RadauCoeff5 &rcoeff5,
            RadauCoeff7 &rcoeff7) {
  /* --- */
  rcoeff3.t311 = .09123239487089294279155;
  rcoeff3.t312 = -.141255295020954208428;
  rcoeff3.t313 = -.03002919410514742449186;
  rcoeff3.t321 = .2417179327071070189575;
  rcoeff3.t322 = .204129352293799931996;
  rcoeff3.t323 = .3829421127572619377954;
  rcoeff3.t331 = .9660481826150929361906;
  rcoeff3.ti311 = 4.325579890063155351024;
  rcoeff3.ti312 = .3391992518158098695428;
  rcoeff3.ti313 = .5417705399358748711865;
  rcoeff3.ti321 = -4.178718591551904727346;
  rcoeff3.ti322 = -.3276828207610623870825;
  rcoeff3.ti323 = .4766235545005504519601;
  rcoeff3.ti331 = -.5028726349457868759512;
  rcoeff3.ti332 = 2.571926949855605429187;
  rcoeff3.ti333 = -.5960392048282249249688;
  if (nsmax <= 3) {
    return;
  }
  rcoeff5.t511 = -.01251758622050104589014;
  rcoeff5.t512 = -.01024204781790882707009;
  rcoeff5.t513 = .04767387729029572386318;
  rcoeff5.t514 = -.01147851525522951470794;
  rcoeff5.t515 = -.01401985889287541028108;
  rcoeff5.t521 = -.001491670151895382429004;
  rcoeff5.t522 = .05017286451737105816299;
  rcoeff5.t523 = -.09433181918161143698066;
  rcoeff5.t524 = -.007668830749180162885157;
  rcoeff5.t525 = .02470857842651852681253;
  rcoeff5.t531 = -.07298187638808714862266;
  rcoeff5.t532 = -.2305395340434179467214;
  rcoeff5.t533 = .1027030453801258997922;
  rcoeff5.t534 = .01939846399882895091122;
  rcoeff5.t535 = .08180035370375117083639;
  rcoeff5.t541 = -.3800914400035681041264;
  rcoeff5.t542 = .3778939022488612495439;
  rcoeff5.t543 = .4667441303324943592896;
  rcoeff5.t544 = .4076011712801990666217;
  rcoeff5.t545 = .1996824278868025259365;
  rcoeff5.t551 = -.9219789736812104884883;
  rcoeff5.ti511 = -30.04156772154440162771;
  rcoeff5.ti512 = -13.86510785627141316518;
  rcoeff5.ti513 = -3.480002774795185561828;
  rcoeff5.ti514 = 1.032008797825263422771;
  rcoeff5.ti515 = -.8043030450739899174753;
  rcoeff5.ti521 = 5.344186437834911598895;
  rcoeff5.ti522 = 4.593615567759161004454;
  rcoeff5.ti523 = -3.036360323459424298646;
  rcoeff5.ti524 = 1.05066019023145886386;
  rcoeff5.ti525 = -.2727786118642962705386;
  rcoeff5.ti531 = 3.748059807439804860051;
  rcoeff5.ti532 = -3.984965736343884667252;
  rcoeff5.ti533 = -1.044415641608018792942;
  rcoeff5.ti534 = 1.184098568137948487231;
  rcoeff5.ti535 = -.4499177701567803688988;
  rcoeff5.ti541 = -33.04188021351900000806;
  rcoeff5.ti542 = -17.37695347906356701945;
  rcoeff5.ti543 = -.1721290632540055611515;
  rcoeff5.ti544 = -.09916977798254264258817;
  rcoeff5.ti545 = .5312281158383066671849;
  rcoeff5.ti551 = -8.6114439798752919777;
  rcoeff5.ti552 = 9.699991409528808231336;
  rcoeff5.ti553 = 1.914728639696874284851;
  rcoeff5.ti554 = 2.418692006084940026427;
  rcoeff5.ti555 = -1.047463487935337418694;
  if (nsmax <= 5) {
    return;
  }
  rcoeff7.t711 = -.002153754627310526422828;
  rcoeff7.t712 = .02156755135132077338691;
  rcoeff7.t713 = .008783567925144144407326;
  rcoeff7.t714 = -.004055161452331023898198;
  rcoeff7.t715 = .004427232753268285479678;
  rcoeff7.t716 = -.001238646187952874056377;
  rcoeff7.t717 = -.002760617480543852499548;
  rcoeff7.t721 = .001600025077880428526831;
  rcoeff7.t722 = -.03813164813441154669442;
  rcoeff7.t723 = -.02152556059400687552385;
  rcoeff7.t724 = .008415568276559589237177;
  rcoeff7.t725 = -.004031949570224549492304;
  rcoeff7.t726 = -6.666635339396338181761e-5;
  rcoeff7.t727 = .003185474825166209848748;
  rcoeff7.t731 = -.00405910730194768309165;
  rcoeff7.t732 = .05739650893938171539757;
  rcoeff7.t733 = .05885052920842679105612;
  rcoeff7.t734 = -.008560431061603432060177;
  rcoeff7.t735 = -.006923212665023908924141;
  rcoeff7.t736 = -.002352180982943338340535;
  rcoeff7.t737 = 4.169077725297562691409e-4;
  rcoeff7.t741 = -.01575048807937684420346;
  rcoeff7.t742 = -.03821469359696835048464;
  rcoeff7.t743 = -.1657368112729438512412;
  rcoeff7.t744 = -.03737124230238445741907;
  rcoeff7.t745 = .008239007298507719404499;
  rcoeff7.t746 = .003115071152346175252726;
  rcoeff7.t747 = .02511660491343882192836;
  rcoeff7.t751 = -.1129776610242208076086;
  rcoeff7.t752 = -.2491742124652636863308;
  rcoeff7.t753 = .2735633057986623212132;
  rcoeff7.t754 = .005366761379181770094279;
  rcoeff7.t755 = .1932111161012620144312;
  rcoeff7.t756 = .1017177324817151468081;
  rcoeff7.t757 = .09504502035604622821039;
  rcoeff7.t761 = -.4583810431839315010281;
  rcoeff7.t762 = .5315846490836284292051;
  rcoeff7.t763 = .4863228366175728940567;
  rcoeff7.t764 = .5265742264584492629141;
  rcoeff7.t765 = .2755343949896258141929;
  rcoeff7.t766 = .5217519452747652852946;
  rcoeff7.t767 = .1280719446355438944141;
  rcoeff7.t771 = -.8813915783538183763135;
  rcoeff7.ti711 = -258.1319263199822292761;
  rcoeff7.ti712 = -189.073763081398508952;
  rcoeff7.ti713 = -49.08731481793013119445;
  rcoeff7.ti714 = -4.110647469661428418112;
  rcoeff7.ti715 = -4.053447889315563304175;
  rcoeff7.ti716 = 3.112755366607346076554;
  rcoeff7.ti717 = -1.646774913558444650169;
  rcoeff7.ti721 = -3.007390169451292131731;
  rcoeff7.ti722 = -11.01586607876577132911;
  rcoeff7.ti723 = 1.487799456131656281486;
  rcoeff7.ti724 = 2.130388159559282459432;
  rcoeff7.ti725 = -1.816141086817565624822;
  rcoeff7.ti726 = 1.134325587895161100083;
  rcoeff7.ti727 = -.414699045943303531993;
  rcoeff7.ti731 = -8.441963188321084681757;
  rcoeff7.ti732 = -.6505252740575150028169;
  rcoeff7.ti733 = 6.940670730369876478804;
  rcoeff7.ti734 = -3.205047525597898431565;
  rcoeff7.ti735 = 1.071280943546478589783;
  rcoeff7.ti736 = -.354850749121622187973;
  rcoeff7.ti737 = .09198549132786554154409;
  rcoeff7.ti741 = 74.67833223502269977153;
  rcoeff7.ti742 = 87.40858897990081640204;
  rcoeff7.ti743 = 4.024158737379997877014;
  rcoeff7.ti744 = -3.714806315158364186639;
  rcoeff7.ti745 = -3.430093985982317350741;
  rcoeff7.ti746 = 2.696604809765312378853;
  rcoeff7.ti747 = -.9386927436075461933568;
  rcoeff7.ti751 = 58.35652885190657724237;
  rcoeff7.ti752 = -10.06877395780018096325;
  rcoeff7.ti753 = -30.36638884256667120811;
  rcoeff7.ti754 = -1.020020865184865985027;
  rcoeff7.ti755 = -.1124175003784249621267;
  rcoeff7.ti756 = 1.8906408310003776228;
  rcoeff7.ti757 = -.9716486393831482282172;
  rcoeff7.ti761 = -299.1862480282520966786;
  rcoeff7.ti762 = -243.0407453687447911819;
  rcoeff7.ti763 = -48.77710407803786921219;
  rcoeff7.ti764 = -2.03867190574193440528;
  rcoeff7.ti765 = 1.673560239861084944268;
  rcoeff7.ti766 = -1.087374032057106164456;
  rcoeff7.ti767 = .9019382492960993738427;
  rcoeff7.ti771 = -93.07650289743530591157;
  rcoeff7.ti772 = 23.88163105628114427703;
  rcoeff7.ti773 = 39.2788807308138438271;
  rcoeff7.ti774 = 14.38891568549108006988;
  rcoeff7.ti775 = -3.510438399399361221087;
  rcoeff7.ti776 = 4.863284885566180701215;
  rcoeff7.ti777 = -2.2464827295912399164;
  return;
}

void coercv(const int ns, double cs[], double dd[], double &u1, double alph[],
            double beta[]) {
  /* System generated locals */
  double d1, d2;

  /* Local variables */
  double sq6, st9, bet, alp, cno;

  /* Function Body */
  cs[0] = 0.;
  cs[ns] = 1.;
  switch (ns) {
  case 1:
    goto L1;
  case 2:
    goto L11;
  case 3:
    goto L3;
  case 4:
    goto L11;
  case 5:
    goto L5;
  case 6:
    goto L11;
  case 7:
    goto L7;
  default:;
  }
L11:
  return;
L1:
  cs[1] = 1.;
  u1 = 1.;
  dd[0] = -1.;
  return;
L3:
  sq6 = sqrt(6.);
  cs[1] = (4. - sq6) / 10.;
  cs[2] = (sq6 + 4.) / 10.;
  st9 = pow(c_b140, c_b141);
  u1 = (st9 * (st9 - 1) + 6.) / 30.;
  alp = (12. - st9 * (st9 - 1)) / 60.;
  bet = st9 * (st9 + 1) * sqrt(3.) / 60.;
  /* Computing 2nd power */
  d1 = alp;
  /* Computing 2nd power */
  d2 = bet;
  cno = d1 * d1 + d2 * d2;
  u1 = 1.0 / u1;
  alph[0] = alp / cno;
  beta[0] = bet / cno;
  return;
L5:
  cs[1] = .05710419611451768219312;
  cs[2] = .27684301363812382768;
  cs[3] = .5835904323689168200567;
  cs[4] = .8602401356562194478479;
  dd[1 - 1] = -27.78093394406463730479;
  dd[2 - 1] = 3.641478498049213152712;
  dd[3 - 1] = -1.252547721169118720491;
  dd[4 - 1] = .5920031671845428725662;
  dd[5 - 1] = -.2;
  u1 = 6.286704751729276645173;
  alph[0] = 3.655694325463572258243;
  beta[0] = 6.543736899360077294021;
  alph[1] = 5.70095329867178941917;
  beta[1] = 3.210265600308549888425;
  return;
L7:
  cs[1] = .02931642715978489197205;
  cs[2] = .14807859966848429185;
  cs[3] = .3369846902811542990971;
  cs[4] = .5586715187715501320814;
  cs[5] = .7692338620300545009169;
  cs[6] = .9269456713197411148519;
  dd[1 - 1] = -54.37443689412861451458;
  dd[2 - 1] = 7.000024004259186512041;
  dd[3 - 1] = -2.355661091987557192256;
  dd[4 - 1] = 1.132289066106134386384;
  dd[5 - 1] = -.6468913267673587118673;
  dd[6 - 1] = .3875333853753523774248;
  dd[7 - 1] = -.1428571428571428571429;
  u1 = 8.936832788405216337302;
  alph[1 - 1] = 4.378693561506806002523;
  beta[1 - 1] = 10.16969328379501162732;
  alph[2 - 1] = 7.141055219187640105775;
  beta[2 - 1] = 6.623045922639275970621;
  alph[3 - 1] = 8.511834825102945723051;
  beta[3 - 1] = 3.281013624325058830036;
}

double contra(const int i, const double x, const Vector<double> &cont,
              const int lrc, const RadauWeight &weight) {

  /* ---------------------------------------------------------- */
  /*     THIS FUNCTION CAN BE USED FOR CONINUOUS OUTPUT. IT PROVIDES AN */
  /*     APPROXIMATION TO THE I-TH COMPONENT OF THE SOLUTION AT X. */
  /*     IT GIVES THE VALUE OF THE COLLOCATION POLYNOMIAL, DEFINED FOR */
  /*     THE LAST SUCCESSFULLY COMPUTED STEP (BY RADAU). */
  /* ---------------------------------------------------------- */

  /* Function Body */
  double s = (x - weight.xsol) / weight.hsol + 1.;
  double ret_val = cont[i + weight.ns * weight.nn - 1];
  for (int k = weight.ns - 1; k >= 0; --k) {
    ret_val = cont[i + k * weight.nn - 1] +
              (s - weight.c[weight.ns - k - 1]) * ret_val;
  }
  return ret_val;
}

int radcov(const size_t n, OdeFun &fcn, double &x, Vector<double> &y,
           const double xend, const double hmax, double h, Vector<double> &rtol,
           Vector<double> &atol, const int ns, OdeJac &jac, S_fp solout,
           bool iout, int &idid, const int nmax, const double uround,
           const double safe, const double thet, const double quot1,
           const double quot2, const int nit1, const int ijob, bool startn,
           const int nind1, const int nind2, const int nind3, const bool pred,
           const double facl, const double facr, const int nsmin,
           const int nsmax, const int nnms, const int nm1ns, const int nmee,
           Vector<double> &zz, Vector<double> &y0, Vector<double> &scal,
           Vector<double> &ff, Matrix<double> &fjac, Matrix<double> &e1,
           Matrix<double> &ee2, Vector<double> &cont, Vector<int> &ip1,
           Vector<int> &ip2, Vector<int>&iphes, double *vitu, double *vitd, double *hhou,
           double *hhod, int nfcn, int njac, int nstep, int naccpt, int nrejct,
           int ndec, int nsol, double *rpar, int *ipar) {
  /* Format strings */
  char fmt_979[] = "(\002 EXIT OF RADAU AT X=\002,e18.4)";

  /* System generated locals */
  int fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e1_dim1, e1_offset,
      ee2_dim1, ee2_offset, i__1, i__2, i__3, i__4;
  double d__1, d__2, d__3, d__4, d__5;

  /* Builtin functions */
  double sqrt(double), d_sign(double *, double *), pow_dd(double *, double *),
      pow_di(double *, int *);
  int s_wsfe(cilist *), do_fio(int *, char *, ftnlen), e_wsfe(),
      s_wsle(cilist *), do_lio(int *, int *, char *, ftnlen), e_wsle();

  /* Local variables */
  int i__, j, k, l;
  double a1, a2, a3;
  int j1, n2, n3, n4, n5, n6;
  double u1, c31, dd[7], c32, ak;
  int md, in, mm;
  double qt, dd1, dd2, dd3, ak1, ak2, ak3, f1i, f2i, f3i, c1q, c2q, c3q, sq6;
  int nns, lrc;
  double fac1;
  int ier, iad;
  double xph, z1i, z2i, c31m1, c32m1, z3i;
  int nit;
  double thq;
  int in2;
  double err, ccq, val, z4i, z5i, z6i, beta[7], alph[7];
  bool last;
  double expo, c31mc2, hold, hopt, xold, rtol1, atol1;
  int lbeg, lend;
  double delt;
  int newt;
  double hhfac, dyno, dyth;
  int ichan, iadd;
  double fact, betan[7], z7i, fac, quot;
  int ikeep;
  double hnew, hacc, alphn[7], denom, theta, ysafe, hmaxn;
  int nsing;
  double expmi, fnewt;
  bool first;
  int nsnew;
  bool unexn, unexp;
  int irtrn, nrsol, nsolu;
  double hquot, xosol, qnewt, quott, acont3;
  bool index1, index2, index3, caljac, change;
  double faccon;
  bool calhes;
  double erracc;
  bool variab;
  int mujacj;
  double facgus;
  bool reject;
  int mujacp;
  double thetat, dynold, posneg;
  double thqold;
  double expmns;

  /* ---------------------------------------------------------- */
  /*     CORE INTEGRATOR FOR RADAU */
  /*     PARAMETERS SAME AS IN RADAU WITH WORKSPACE ADDED */
  /* ---------------------------------------------------------- */
  /*         DECLARATIONS */
  /* ---------------------------------------------------------- */
  /* --- THIS PARAMETER HAS TO BE CHANGED IF NUMBER OF STAGES IS >=7 */
  /* *** *** *** *** *** *** *** */
  /*  INITIALISATIONS */
  /* *** *** *** *** *** *** *** */
  /* -------- CHECK THE INDEX OF THE PROBLEM ----- */
  /* ------- COMPUTE MASS MATRIX FOR IMPLICIT CASE ---------- */
  variab = nsmin < nsmax;
  /* ---------- CONSTANTS --------- */
  expo = 1. / (ns + 1.);
  sq6 = sqrt(6.);
  c31 = (4. - sq6) / 10.;
  c32 = (sq6 + 4.) / 10.;
  c31m1 = c31 - 1.;
  c32m1 = c32 - 1.;
  c31mc2 = c31 - c32;
  dd1 = -(sq6 * 7. + 13.) / 3.;
  dd2 = (sq6 * 7. - 13.) / 3.;
  dd3 = -.33333333333333331;
  n2 = n << 1;
  n3 = n * 3;
  n4 = n << 2;
  n5 = n * 5;
  n6 = n * 6;
  unexp = false;
  unexn = false;
  change = false;
  ikeep = 0;
  ichan = 0;
  theta = 0.;
  thetat = 0.;

  RadauWeight weight{};
  RadauCoeff3 coeff3{};
  RadauCoeff5 coeff5{};
  RadauCoeff7 coeff7{};

  weight.nn = n;
  nns = n * ns;
  weight.ns = ns;
  lrc = nns + n;
  coertv(nsmax, coeff3, coeff5, coeff7);
  coercv(ns, weight.c, dd, u1, alph, beta);
  d__1 = xend - x;
  posneg = std::copysign(c_b90, d__1);
  /* Computing MIN */
  d__2 = abs(hmax), d__3 = (d__1 = xend - x, abs(d__1));
  hmaxn = std::min(d__2, d__3);
  if (abs(h) <= uround * 10.) {
    h = 1e-6;
  }
  /* Computing MIN */
  d__1 = abs(h);
  h = std::min(d__1, hmaxn);
  h = std::copysign(h, posneg);
  hold = h;
  reject = false;
  first = true;
  last = false;
  if ((x + h * 1.0001 - xend) * posneg >= 0.) {
    h = xend - x;
    last = true;
  }
  hopt = h;
  faccon = 1.;
  nsing = 0;
  xold = x;
  if (iout != 0) {
    irtrn = 1;
    nrsol = 1;
    xosol = xold;
    weight.xsol = x;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      cont[i__] = y[i__];
    }
    nsolu = n;
    weight.hsol = hold;
    (*solout)(&nrsol, &xosol, &weight.xsol, &y[1], &cont[1], &lrc, &nsolu,
              &rpar[1], &ipar[1], &irtrn);
    if (irtrn < 0) {
      goto L179;
    }
  }
  expmns = (ns + 1.) / (ns * 2.);
  quott = atol[1] / rtol[1];
  i__1 = n;
  for (i__ = 1; i__ <= i__1; ++i__) {
    quott = atol[i__] / rtol[i__];
    rtol1 = std::pow(rtol[i__], expmns) * 0.1;
    atol1 = rtol1 * quott;
    scal[i__] = atol1 + rtol1 * (d__1 = y[i__], abs(d__1));
  }
  hhfac = h;
  fcn(x, y, y0);
  ++(nfcn);
/* --- BASIC INTEGRATION STEP */
L10:
  /* *** *** *** *** *** *** *** */
  /*  COMPUTATION OF THE JACOBIAN */
  /* *** *** *** *** *** *** *** */
  ++(njac);
  /* --- COMPUTE JACOBIAN MATRIX ANALYTICALLY */
  jac(x, y, fjac);
  caljac = true;
  calhes = true;
L20:
  /* --- CHANGE THE ORDER HERE IF NECESSARY */
  if (variab) {
    nsnew = ns;
    ++ichan;
    hquot = h / hold;
    /* Computing MIN */
    /* Computing MAX */
    d__3 = theta, d__4 = thetat * 0.5;
    d__1 = 10.0, d__2 = std::max(d__3, d__4);
    thetat = std::min(d__1, d__2);
    if (newt > 1 && thetat <= *vitu && hquot < *hhou && hquot > *hhod) {
      /* Computing MIN */
      i__1 = nsmax, i__2 = ns + 2;
      nsnew = std::min(i__1, i__2);
    }
    if (thetat >= *vitd || unexp) {
      /* Computing MAX */
      i__1 = nsmin, i__2 = ns - 2;
      nsnew = std::max(i__1, i__2);
    }
    if (ichan >= 1 && unexn) {
      /* Computing MAX */
      i__1 = nsmin, i__2 = ns - 2;
      nsnew = std::max(i__1, i__2);
    }
    if (ichan <= 10) {
      nsnew = std::min(ns, nsnew);
    }
    change = ns != nsnew;
    unexn = false;
    unexp = false;
    if (change) {
      ns = nsnew;
      ichan = 1;
      nns = n * ns;
      weight.ns = ns;
      lrc = nns + n;
      coercv(ns, weight.c, dd, u1, alph, beta);
      expo = 1.0 / (ns + 1.0);
      expmns = (ns + 1.0) / (ns * 2.0);
      rtol1 = pow(rtol[1], expmns) * 0.1;
      atol1 = rtol1 * quott;
      i__1 = n;
      for (i__ = 1; i__ <= i__1; ++i__) {
        quott = atol[i__] / rtol[i__];
        rtol1 = pow_dd(&rtol[i__], &expmns) * .1;
        atol1 = rtol1 * quott;
        scal[i__] = atol1 + rtol1 * (d__1 = y[i__], abs(d__1));
      }
    }
  }
  /* --- COMPUTE THE MATRICES E1 AND E2 AND THEIR DECOMPOSITIONS */
  fac1 = u1 / h;
  ier = decomr(fjac, fac1, e1, ip1);
  if (ier != 0) {
    goto L78;
  }
  i__1 = (ns - 1) / 2;
  for (k = 1; k <= i__1; ++k) {
    alphn[k - 1] = alph[k - 1] / h;
    betan[k - 1] = beta[k - 1] / h;
    iad = ((k - 1) << 1) * nm1 + 1;
    ier = decomc(fjac, alphn[k - 1], betan[k - 1],
            ee2[iad * ee2_dim1 + 1], ee2[(iad + nm1) * ee2_dim1 + 1],
            ip2[(k - 1) * nm1 + 1]);
    if (ier != 0) {
      goto L78;
    }
  }
  ++(ndec);
L30:
  if (variab && ikeep == 1) {
    ++ichan;
    ikeep = 0;
    if (ichan >= 10 && ns < nsmax) {
      goto L20;
    }
  }
  ++(nstep);
  if (nstep > nmax) {
    goto L178;
  }
  if (abs(h) * 0.1 <= abs(x) * uround) {
    goto L177;
  }
  if (index2) {
    i__1 = nind1 + nind2;
    for (i__ = nind1 + 1; i__ <= i__1; ++i__) {
      scal[i__] /= hhfac;
    }
  }
  if (index3) {
    i__1 = nind1 + nind2 + nind3;
    for (i__ = nind1 + nind2 + 1; i__ <= i__1; ++i__) {
      scal[i__] /= hhfac * hhfac;
    }
  }
  xph = x + h;
  /* *** *** *** *** *** *** *** */
  /* *** *** *** *** *** *** *** */
  if (ns == 3) {
    /* *** *** *** *** *** *** *** */
    /* *** *** *** *** *** *** *** */
    if (first || startn || change) {
      i__1 = nns;
      for (i__ = 1; i__ <= i__1; ++i__) {
        zz[i__] = 0.;
        ff[i__] = 0.;
      }
    } else {
      hquot = h / hold;
      c3q = hquot;
      c1q = c31 * c3q;
      c2q = c32 * c3q;
      i__1 = n;
      for (i__ = 1; i__ <= i__1; ++i__) {
        ak1 = cont[i__ + n];
        ak2 = cont[i__ + n2];
        ak3 = cont[i__ + n3];
        z1i = c1q * (ak1 + (c1q - c32m1) * (ak2 + (c1q - c31m1) * ak3));
        z2i = c2q * (ak1 + (c2q - c32m1) * (ak2 + (c2q - c31m1) * ak3));
        z3i = c3q * (ak1 + (c3q - c32m1) * (ak2 + (c3q - c31m1) * ak3));
        zz[i__] = z1i;
        zz[i__ + n] = z2i;
        zz[i__ + n2] = z3i;
        ff[i__] =
            coeff3.ti311 * z1i + coeff3.ti312 * z2i + coeff3.ti313 * z3i;
        ff[i__ + n] =
            coeff3.ti321 * z1i + coeff3.ti322 * z2i + coeff3.ti323 * z3i;
        ff[i__ + n2] =
            coeff3.ti331 * z1i + coeff3.ti332 * z2i + coeff3.ti333 * z3i;
      }
    }
    /* *** *** *** *** *** *** *** */
    /*  LOOP FOR THE SIMPLIFIED NEWTON ITERATION */
    /* *** *** *** *** *** *** *** */
    newt = 0;
    nit = nit1;
    expmi = 1. / expmns;
    /* Computing MAX */
    /* Computing MIN */
    d__5 = expmi - 1.;
    d__3 = .03, d__4 = pow_dd(&rtol1, &d__5);
    d__1 = uround * 10 / rtol1, d__2 = std::min(d__3, d__4);
    fnewt = std::max(d__1, d__2);
    d__1 = std::max(faccon, uround);
    faccon = std::pow(d__1, c_b100);
    theta = std::abs(thet);
  L40:
    if (newt >= nit) {
      goto L78;
    }
    /* ---     COMPUTE THE RIGHT-HAND SIDE */
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      cont[i__] = y[i__] + zz[i__];
    }
    d__1 = x + c31 * h;
    fcn(d__1, cont, zz);
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      cont[i__] = y[i__] + zz[i__ + n];
    }
    d__1 = x + c32 * h;
    fcn(d__1, cont, zz[n + 1], &rpar[1], &ipar[1]);
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      cont[i__] = y[i__] + zz[i__ + n2];
    }
    (*fcn)(n, &xph, &cont[1], &zz[n2 + 1], &rpar[1], &ipar[1]);
    nfcn += 3;
    /* ---     SOLVE THE LINEAR SYSTEMS */
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      a1 = zz[i__];
      a2 = zz[i__ + n];
      a3 = zz[i__ + n2];
      zz[i__] = rcoeff3.ti311 * a1 + rcoeff3.ti312 * a2 + rcoeff3.ti313 * a3;
      zz[i__ + n] =
          rcoeff3.ti321 * a1 + rcoeff3.ti322 * a2 + rcoeff3.ti323 * a3;
      zz[i__ + n2] =
          rcoeff3.ti331 * a1 + rcoeff3.ti332 * a2 + rcoeff3.ti333 * a3;
    }
    slvrad_(n, &fjac[fjac_offset], ldjac, mljac, mujac, &fmas[fmas_offset],
            ldmas, mlmas, mumas, m1, m2, nm1, &fac1, alphn, betan,
            &e1[e1_offset], &ee2[ee2_offset], &ee2[(nm1 + 1) * ee2_dim1 + 1],
            lde1, &zz[1], &zz[n + 1], &zz[n2 + 1], &ff[1], &ff[n + 1],
            &ff[n2 + 1], &cont[1], &ip1[1], &ip2[1], &iphes[1], &ier, ijob);
    ++(nsol);
    ++newt;
    dyno = 0.;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      denom = scal[i__];
      /* Computing 2nd power */
      d__1 = zz[i__] / denom;
      /* Computing 2nd power */
      d__2 = zz[i__ + n] / denom;
      /* Computing 2nd power */
      d__3 = zz[i__ + n2] / denom;
      dyno = dyno + d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
    }
    dyno = sqrt(dyno / nns);
    /* ---     BAD CONVERGENCE OR NUMBER OF ITERATIONS TO LARGE */
    if (newt > 1 && newt < nit) {
      thq = dyno / dynold;
      if (newt == 2) {
        theta = thq;
      } else {
        theta = sqrt(thq * thqold);
      }
      thqold = thq;
      if (theta < .99) {
        faccon = theta / (1. - theta);
        i__1 = nit - 1 - newt;
        dyth = faccon * dyno * pow_di(&theta, &i__1) / fnewt;
        if (dyth >= 1.) {
          /* Computing MAX */
          d__1 = 1e-4, d__2 = min(20., dyth);
          qnewt = max(d__1, d__2);
          d__1 = -1. / (nit + 4. - 1 - newt);
          hhfac = pow_dd(&qnewt, &d__1) * .8;
          *h = hhfac * *h;
          reject = true;
          last = false;
          if (hhfac <= .5) {
            unexn = true;
          }
          if (caljac) {
            goto L20;
          }
          goto L10;
        }
      } else {
        goto L78;
      }
    }
    dynold = max(dyno, *uround);
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      in = i__ + n;
      in2 = i__ + n2;
      f1i = ff[i__] + zz[i__];
      f2i = ff[in] + zz[in];
      f3i = ff[in2] + zz[in2];
      ff[i__] = f1i;
      ff[in] = f2i;
      ff[in2] = f3i;
      zz[i__] = rcoeff3.t311 * f1i + rcoeff3.t312 * f2i + rcoeff3.t313 * f3i;
      zz[in] = rcoeff3.t321 * f1i + rcoeff3.t322 * f2i + rcoeff3.t323 * f3i;
      zz[in2] = rcoeff3.t331 * f1i + f2i;
    }
    if (faccon * dyno > fnewt) {
      goto L40;
    }
    /* --- ERROR ESTIMATION */
    estrad_(n, &fjac[fjac_offset], ldjac, mljac, mujac, &fmas[fmas_offset],
            ldmas, mlmas, mumas, h, &dd1, &dd2, &dd3, (S_fp)fcn, nfcn, &y0[1],
            &y[1], ijob, x, m1, m2, nm1, &e1[e1_offset], lde1, &zz[1],
            &zz[n + 1], &zz[n2 + 1], &cont[1], &ff[1], &ff[n + 1], &ip1[1],
            &iphes[1], &scal[1], &err, &first, &reject, &fac1, &rpar[1],
            &ipar[1]);
    /*       --- COMPUTE FINITE DIFFERENCES FOR DENSE OUTPUT */
    if (err < 1.) {
      i__1 = n;
      for (i__ = 1; i__ <= i__1; ++i__) {
        y[i__] += zz[i__ + n2];
        z2i = zz[i__ + n];
        z1i = zz[i__];
        cont[i__ + n] = (z2i - zz[i__ + n2]) / c32m1;
        ak = (z1i - z2i) / c31mc2;
        acont3 = z1i / c31;
        acont3 = (ak - acont3) / c32;
        cont[i__ + n2] = (ak - cont[i__ + n]) / c31m1;
        cont[i__ + n3] = cont[i__ + n2] - acont3;
      }
    }
    /* *** *** *** *** *** *** *** */
    /* *** *** *** *** *** *** *** */
  } else {
    if (ns == 5) {
      /* *** *** *** *** *** *** *** */
      /* *** *** *** *** *** *** *** */
      if (first || *startn || change) {
        i__1 = nns;
        for (i__ = 1; i__ <= i__1; ++i__) {
          zz[i__] = 0.;
          ff[i__] = 0.;
        }
      } else {
        hquot = *h / hold;
        i__1 = ns;
        for (k = 1; k <= i__1; ++k) {
          ccq = weight_1.c__[k] * hquot;
          i__2 = n;
          for (i__ = 1; i__ <= i__2; ++i__) {
            val = cont[i__ + ns * n];
            for (l = ns - 1; l >= 1; --l) {
              val = cont[i__ + l * n] + (ccq - weight_1.c__[ns - l] + 1.) * val;
            }
            zz[i__ + (k - 1) * n] = ccq * val;
          }
        }
        i__1 = n;
        for (i__ = 1; i__ <= i__1; ++i__) {
          z1i = zz[i__];
          z2i = zz[i__ + n];
          z3i = zz[i__ + n2];
          z4i = zz[i__ + n3];
          z5i = zz[i__ + n4];
          ff[i__] = rcoeff5.ti511 * z1i + rcoeff5.ti512 * z2i +
                    rcoeff5.ti513 * z3i + rcoeff5.ti514 * z4i +
                    rcoeff5.ti515 * z5i;
          ff[i__ + n] = rcoeff5.ti521 * z1i + rcoeff5.ti522 * z2i +
                        rcoeff5.ti523 * z3i + rcoeff5.ti524 * z4i +
                        rcoeff5.ti525 * z5i;
          ff[i__ + n2] = rcoeff5.ti531 * z1i + rcoeff5.ti532 * z2i +
                         rcoeff5.ti533 * z3i + rcoeff5.ti534 * z4i +
                         rcoeff5.ti535 * z5i;
          ff[i__ + n3] = rcoeff5.ti541 * z1i + rcoeff5.ti542 * z2i +
                         rcoeff5.ti543 * z3i + rcoeff5.ti544 * z4i +
                         rcoeff5.ti545 * z5i;
          ff[i__ + n4] = rcoeff5.ti551 * z1i + rcoeff5.ti552 * z2i +
                         rcoeff5.ti553 * z3i + rcoeff5.ti554 * z4i +
                         rcoeff5.ti555 * z5i;
        }
      }
      /* *** *** *** *** *** *** *** */
      /*  LOOP FOR THE SIMPLIFIED NEWTON ITERATION */
      /* *** *** *** *** *** *** *** */
      newt = 0;
      nit = nit1 + 5;
      expmi = 1. / expmns;
      /* Computing MAX */
      /* Computing MIN */
      d__5 = expmi - 1.;
      d__3 = .03, d__4 = pow_dd(&rtol1, &d__5);
      d__1 = *uround * 10 / rtol1, d__2 = min(d__3, d__4);
      fnewt = max(d__1, d__2);
      d__1 = max(faccon, *uround);
      faccon = pow_dd(&d__1, &c_b100);
      theta = abs(*thet);
    L140:
      if (newt >= nit) {
        goto L78;
      }
      /* ---     COMPUTE THE RIGHT-HAND SIDE */
      i__1 = ns - 1;
      for (k = 0; k <= i__1; ++k) {
        iadd = k * n;
        i__2 = n;
        for (i__ = 1; i__ <= i__2; ++i__) {
          cont[i__] = y[i__] + zz[iadd + i__];
        }
        d__1 = *x + weight_1.c__[k + 1] * *h;
        (*fcn)(n, &d__1, &cont[1], &zz[iadd + 1], &rpar[1], &ipar[1]);
      }
      nfcn += ns;
      /* ---     SOLVE THE LINEAR SYSTEMS */
      i__1 = n;
      for (i__ = 1; i__ <= i__1; ++i__) {
        z1i = zz[i__];
        z2i = zz[i__ + n];
        z3i = zz[i__ + n2];
        z4i = zz[i__ + n3];
        z5i = zz[i__ + n4];
        zz[i__] = rcoeff5.ti511 * z1i + rcoeff5.ti512 * z2i +
                  rcoeff5.ti513 * z3i + rcoeff5.ti514 * z4i +
                  rcoeff5.ti515 * z5i;
        zz[i__ + n] = rcoeff5.ti521 * z1i + rcoeff5.ti522 * z2i +
                      rcoeff5.ti523 * z3i + rcoeff5.ti524 * z4i +
                      rcoeff5.ti525 * z5i;
        zz[i__ + n2] = rcoeff5.ti531 * z1i + rcoeff5.ti532 * z2i +
                       rcoeff5.ti533 * z3i + rcoeff5.ti534 * z4i +
                       rcoeff5.ti535 * z5i;
        zz[i__ + n3] = rcoeff5.ti541 * z1i + rcoeff5.ti542 * z2i +
                       rcoeff5.ti543 * z3i + rcoeff5.ti544 * z4i +
                       rcoeff5.ti545 * z5i;
        zz[i__ + n4] = rcoeff5.ti551 * z1i + rcoeff5.ti552 * z2i +
                       rcoeff5.ti553 * z3i + rcoeff5.ti554 * z4i +
                       rcoeff5.ti555 * z5i;
      }
      slvrar_(n, &fjac[fjac_offset], ldjac, mljac, mujac, &fmas[fmas_offset],
              ldmas, mlmas, mumas, m1, m2, nm1, &fac1, &e1[e1_offset], lde1,
              &zz[1], &ff[1], &ip1[1], &iphes[1], &ier, ijob);
      for (k = 1; k <= 2; ++k) {
        iad = (k - 1 << 1) * nm1 + 1;
        slvrai_(n, &fjac[fjac_offset], ldjac, mljac, mujac, &fmas[fmas_offset],
                ldmas, mlmas, mumas, m1, m2, nm1, &alphn[k - 1], &betan[k - 1],
                &ee2[iad * ee2_dim1 + 1], &ee2[(iad + nm1) * ee2_dim1 + 1],
                lde1, &zz[((k << 1) - 1) * n + 1], &zz[(k << 1) * n + 1],
                &ff[((k << 1) - 1) * n + 1], &ff[(k << 1) * n + 1], &cont[1],
                &ip2[(k - 1) * nm1 + 1], &iphes[1], &ier, ijob);
      }
      ++(nsol);
      ++newt;
      dyno = 0.;
      i__1 = n;
      for (i__ = 1; i__ <= i__1; ++i__) {
        denom = scal[i__];
        i__2 = ns - 1;
        for (k = 0; k <= i__2; ++k) {
          /* Computing 2nd power */
          d__1 = zz[i__ + k * n] / denom;
          dyno += d__1 * d__1;
        }
      }
      dyno = sqrt(dyno / nns);
      /* ---     BAD CONVERGENCE OR NUMBER OF ITERATIONS TO LARGE */
      if (newt > 1 && newt < nit) {
        thq = dyno / dynold;
        if (newt == 2) {
          theta = thq;
        } else {
          theta = sqrt(thq * thqold);
        }
        thqold = thq;
        if (theta < .99) {
          faccon = theta / (1. - theta);
          i__1 = nit - 1 - newt;
          dyth = faccon * dyno * pow_di(&theta, &i__1) / fnewt;
          if (dyth >= 1.) {
            /* Computing MAX */
            d__1 = 1e-4, d__2 = min(20., dyth);
            qnewt = max(d__1, d__2);
            d__1 = -1. / (nit + 4. - 1 - newt);
            hhfac = pow_dd(&qnewt, &d__1) * .8;
            *h = hhfac * *h;
            reject = true;
            last = false;
            if (hhfac <= .5) {
              unexn = true;
            }
            if (caljac) {
              goto L20;
            }
            goto L10;
          }
        } else {
          goto L78;
        }
      }
      dynold = max(dyno, *uround);
      i__1 = n;
      for (i__ = 1; i__ <= i__1; ++i__) {
        z1i = ff[i__] + zz[i__];
        z2i = ff[i__ + n] + zz[i__ + n];
        z3i = ff[i__ + n2] + zz[i__ + n2];
        z4i = ff[i__ + n3] + zz[i__ + n3];
        z5i = ff[i__ + n4] + zz[i__ + n4];
        ff[i__] = z1i;
        ff[i__ + n] = z2i;
        ff[i__ + n2] = z3i;
        ff[i__ + n3] = z4i;
        ff[i__ + n4] = z5i;
        zz[i__] = rcoeff5.t511 * z1i + rcoeff5.t512 * z2i + rcoeff5.t513 * z3i +
                  rcoeff5.t514 * z4i + rcoeff5.t515 * z5i;
        zz[i__ + n] = rcoeff5.t521 * z1i + rcoeff5.t522 * z2i +
                      rcoeff5.t523 * z3i + rcoeff5.t524 * z4i +
                      rcoeff5.t525 * z5i;
        zz[i__ + n2] = rcoeff5.t531 * z1i + rcoeff5.t532 * z2i +
                       rcoeff5.t533 * z3i + rcoeff5.t534 * z4i +
                       rcoeff5.t535 * z5i;
        zz[i__ + n3] = rcoeff5.t541 * z1i + rcoeff5.t542 * z2i +
                       rcoeff5.t543 * z3i + rcoeff5.t544 * z4i +
                       rcoeff5.t545 * z5i;
        zz[i__ + n4] = rcoeff5.t551 * z1i + z2i + z4i;
      }
      if (faccon * dyno > fnewt) {
        goto L140;
      }
      /* --- ERROR ESTIMATION */
      estrav_(n, &fjac[fjac_offset], ldjac, mljac, mujac, &fmas[fmas_offset],
              ldmas, mlmas, mumas, h, dd, (S_fp)fcn, nfcn, &y0[1], &y[1], ijob,
              x, m1, m2, nm1, ns, &nns, &e1[e1_offset], lde1, &zz[1], &cont[1],
              &ff[1], &ip1[1], &iphes[1], &scal[1], &err, &first, &reject,
              &fac1, &rpar[1], &ipar[1]);
      /*       --- COMPUTE FINITE DIFFERENCES FOR DENSE OUTPUT */
      if (err < 1.) {
        i__1 = n;
        for (i__ = 1; i__ <= i__1; ++i__) {
          y[i__] += zz[i__ + n4];
          cont[i__ + n5] = zz[i__] / weight_1.c__[1];
        }
        i__1 = ns - 1;
        for (k = 1; k <= i__1; ++k) {
          fact = 1. / (weight_1.c__[ns - k] - weight_1.c__[ns - k + 1]);
          i__2 = n;
          for (i__ = 1; i__ <= i__2; ++i__) {
            cont[i__ + k * n] =
                (zz[i__ + (ns - k - 1) * n] - zz[i__ + (ns - k) * n]) * fact;
          }
        }
        i__1 = ns;
        for (j = 2; j <= i__1; ++j) {
          i__2 = j;
          for (k = ns; k >= i__2; --k) {
            fact = 1. / (weight_1.c__[ns - k] - weight_1.c__[ns - k + j]);
            i__3 = n;
            for (i__ = 1; i__ <= i__3; ++i__) {
              cont[i__ + k * n] =
                  (cont[i__ + k * n] - cont[i__ + (k - 1) * n]) * fact;
            }
          }
        }
      }
      /* *** *** *** *** *** *** *** */
      /* *** *** *** *** *** *** *** */
    } else {
      if (ns == 7) {
        /* *** *** *** *** *** *** *** */
        /* *** *** *** *** *** *** *** */
        if (first || *startn || change) {
          i__1 = nns;
          for (i__ = 1; i__ <= i__1; ++i__) {
            zz[i__] = 0.;
            ff[i__] = 0.;
          }
        } else {
          hquot = *h / hold;
          i__1 = ns;
          for (k = 1; k <= i__1; ++k) {
            ccq = weight_1.c__[k] * hquot;
            i__2 = n;
            for (i__ = 1; i__ <= i__2; ++i__) {
              val = cont[i__ + ns * n];
              for (l = ns - 1; l >= 1; --l) {
                val =
                    cont[i__ + l * n] + (ccq - weight_1.c__[ns - l] + 1.) * val;
              }
              zz[i__ + (k - 1) * n] = ccq * val;
            }
          }
          i__1 = n;
          for (i__ = 1; i__ <= i__1; ++i__) {
            z1i = zz[i__];
            z2i = zz[i__ + n];
            z3i = zz[i__ + n2];
            z4i = zz[i__ + n3];
            z5i = zz[i__ + n4];
            z6i = zz[i__ + n5];
            z7i = zz[i__ + n6];
            ff[i__] = rcoeff7.ti711 * z1i + rcoeff7.ti712 * z2i +
                      rcoeff7.ti713 * z3i + rcoeff7.ti714 * z4i +
                      rcoeff7.ti715 * z5i + rcoeff7.ti716 * z6i +
                      rcoeff7.ti717 * z7i;
            ff[i__ + n] = rcoeff7.ti721 * z1i + rcoeff7.ti722 * z2i +
                          rcoeff7.ti723 * z3i + rcoeff7.ti724 * z4i +
                          rcoeff7.ti725 * z5i + rcoeff7.ti726 * z6i +
                          rcoeff7.ti727 * z7i;
            ff[i__ + n2] = rcoeff7.ti731 * z1i + rcoeff7.ti732 * z2i +
                           rcoeff7.ti733 * z3i + rcoeff7.ti734 * z4i +
                           rcoeff7.ti735 * z5i + rcoeff7.ti736 * z6i +
                           rcoeff7.ti737 * z7i;
            ff[i__ + n3] = rcoeff7.ti741 * z1i + rcoeff7.ti742 * z2i +
                           rcoeff7.ti743 * z3i + rcoeff7.ti744 * z4i +
                           rcoeff7.ti745 * z5i + rcoeff7.ti746 * z6i +
                           rcoeff7.ti747 * z7i;
            ff[i__ + n4] = rcoeff7.ti751 * z1i + rcoeff7.ti752 * z2i +
                           rcoeff7.ti753 * z3i + rcoeff7.ti754 * z4i +
                           rcoeff7.ti755 * z5i + rcoeff7.ti756 * z6i +
                           rcoeff7.ti757 * z7i;
            ff[i__ + n5] = rcoeff7.ti761 * z1i + rcoeff7.ti762 * z2i +
                           rcoeff7.ti763 * z3i + rcoeff7.ti764 * z4i +
                           rcoeff7.ti765 * z5i + rcoeff7.ti766 * z6i +
                           rcoeff7.ti767 * z7i;
            ff[i__ + n6] = rcoeff7.ti771 * z1i + rcoeff7.ti772 * z2i +
                           rcoeff7.ti773 * z3i + rcoeff7.ti774 * z4i +
                           rcoeff7.ti775 * z5i + rcoeff7.ti776 * z6i +
                           rcoeff7.ti777 * z7i;
          }
        }
        /* *** *** *** *** *** *** *** */
        /*  LOOP FOR THE SIMPLIFIED NEWTON ITERATION */
        /* *** *** *** *** *** *** *** */
        newt = 0;
        nit = nit1 + 10;
        expmi = 1. / expmns;
        /* Computing MAX */
        /* Computing MIN */
        d__5 = expmi - 1.;
        d__3 = .03, d__4 = pow_dd(&rtol1, &d__5);
        d__1 = *uround * 10 / rtol1, d__2 = min(d__3, d__4);
        fnewt = max(d__1, d__2);
        d__1 = max(faccon, *uround);
        faccon = pow_dd(&d__1, &c_b100);
        theta = abs(*thet);
      L240:
        if (newt >= nit) {
          goto L78;
        }
        /* ---     COMPUTE THE RIGHT-HAND SIDE */
        i__1 = ns - 1;
        for (k = 0; k <= i__1; ++k) {
          iadd = k * n;
          i__2 = n;
          for (i__ = 1; i__ <= i__2; ++i__) {
            cont[i__] = y[i__] + zz[iadd + i__];
          }
          d__1 = *x + weight_1.c__[k + 1] * *h;
          (*fcn)(n, &d__1, &cont[1], &zz[iadd + 1], &rpar[1], &ipar[1]);
        }
        nfcn += ns;
        /* ---     SOLVE THE LINEAR SYSTEMS */
        i__1 = n;
        for (i__ = 1; i__ <= i__1; ++i__) {
          z1i = zz[i__];
          z2i = zz[i__ + n];
          z3i = zz[i__ + n2];
          z4i = zz[i__ + n3];
          z5i = zz[i__ + n4];
          z6i = zz[i__ + n5];
          z7i = zz[i__ + n6];
          zz[i__] = rcoeff7.ti711 * z1i + rcoeff7.ti712 * z2i +
                    rcoeff7.ti713 * z3i + rcoeff7.ti714 * z4i +
                    rcoeff7.ti715 * z5i + rcoeff7.ti716 * z6i +
                    rcoeff7.ti717 * z7i;
          zz[i__ + n] = rcoeff7.ti721 * z1i + rcoeff7.ti722 * z2i +
                        rcoeff7.ti723 * z3i + rcoeff7.ti724 * z4i +
                        rcoeff7.ti725 * z5i + rcoeff7.ti726 * z6i +
                        rcoeff7.ti727 * z7i;
          zz[i__ + n2] = rcoeff7.ti731 * z1i + rcoeff7.ti732 * z2i +
                         rcoeff7.ti733 * z3i + rcoeff7.ti734 * z4i +
                         rcoeff7.ti735 * z5i + rcoeff7.ti736 * z6i +
                         rcoeff7.ti737 * z7i;
          zz[i__ + n3] = rcoeff7.ti741 * z1i + rcoeff7.ti742 * z2i +
                         rcoeff7.ti743 * z3i + rcoeff7.ti744 * z4i +
                         rcoeff7.ti745 * z5i + rcoeff7.ti746 * z6i +
                         rcoeff7.ti747 * z7i;
          zz[i__ + n4] = rcoeff7.ti751 * z1i + rcoeff7.ti752 * z2i +
                         rcoeff7.ti753 * z3i + rcoeff7.ti754 * z4i +
                         rcoeff7.ti755 * z5i + rcoeff7.ti756 * z6i +
                         rcoeff7.ti757 * z7i;
          zz[i__ + n5] = rcoeff7.ti761 * z1i + rcoeff7.ti762 * z2i +
                         rcoeff7.ti763 * z3i + rcoeff7.ti764 * z4i +
                         rcoeff7.ti765 * z5i + rcoeff7.ti766 * z6i +
                         rcoeff7.ti767 * z7i;
          zz[i__ + n6] = rcoeff7.ti771 * z1i + rcoeff7.ti772 * z2i +
                         rcoeff7.ti773 * z3i + rcoeff7.ti774 * z4i +
                         rcoeff7.ti775 * z5i + rcoeff7.ti776 * z6i +
                         rcoeff7.ti777 * z7i;
        }
        slvrar_(n, &fjac[fjac_offset], ldjac, mljac, mujac, &fmas[fmas_offset],
                ldmas, mlmas, mumas, m1, m2, nm1, &fac1, &e1[e1_offset], lde1,
                &zz[1], &ff[1], &ip1[1], &iphes[1], &ier, ijob);
        for (k = 1; k <= 3; ++k) {
          iad = (k - 1 << 1) * nm1 + 1;
          slvrai_(n, &fjac[fjac_offset], ldjac, mljac, mujac,
                  &fmas[fmas_offset], ldmas, mlmas, mumas, m1, m2, nm1,
                  &alphn[k - 1], &betan[k - 1], &ee2[iad * ee2_dim1 + 1],
                  &ee2[(iad + nm1) * ee2_dim1 + 1], lde1,
                  &zz[((k << 1) - 1) * n + 1], &zz[(k << 1) * n + 1],
                  &ff[((k << 1) - 1) * n + 1], &ff[(k << 1) * n + 1], &cont[1],
                  &ip2[(k - 1) * nm1 + 1], &iphes[1], &ier, ijob);
        }
        ++(nsol);
        ++newt;
        dyno = 0.;
        i__1 = n;
        for (i__ = 1; i__ <= i__1; ++i__) {
          denom = scal[i__];
          i__2 = ns - 1;
          for (k = 0; k <= i__2; ++k) {
            /* Computing 2nd power */
            d__1 = zz[i__ + k * n] / denom;
            dyno += d__1 * d__1;
          }
        }
        dyno = sqrt(dyno / nns);
        /* ---     BAD CONVERGENCE OR NUMBER OF ITERATIONS TO LARGE */
        if (newt > 1 && newt < nit) {
          thq = dyno / dynold;
          if (newt == 2) {
            theta = thq;
          } else {
            theta = sqrt(thq * thqold);
          }
          thqold = thq;
          if (theta < .99) {
            faccon = theta / (1. - theta);
            i__1 = nit - 1 - newt;
            dyth = faccon * dyno * pow_di(&theta, &i__1) / fnewt;
            if (dyth >= 1.) {
              /* Computing MAX */
              d__1 = 1e-4, d__2 = min(20., dyth);
              qnewt = max(d__1, d__2);
              d__1 = -1. / (nit + 4. - 1 - newt);
              hhfac = pow_dd(&qnewt, &d__1) * .8;
              *h = hhfac * *h;
              reject = true;
              last = false;
              if (hhfac <= .5) {
                unexn = true;
              }
              if (caljac) {
                goto L20;
              }
              goto L10;
            }
          } else {
            goto L78;
          }
        }
        dynold = max(dyno, *uround);
        i__1 = n;
        for (i__ = 1; i__ <= i__1; ++i__) {
          z1i = ff[i__] + zz[i__];
          z2i = ff[i__ + n] + zz[i__ + n];
          z3i = ff[i__ + n2] + zz[i__ + n2];
          z4i = ff[i__ + n3] + zz[i__ + n3];
          z5i = ff[i__ + n4] + zz[i__ + n4];
          z6i = ff[i__ + n5] + zz[i__ + n5];
          z7i = ff[i__ + n6] + zz[i__ + n6];
          ff[i__] = z1i;
          ff[i__ + n] = z2i;
          ff[i__ + n2] = z3i;
          ff[i__ + n3] = z4i;
          ff[i__ + n4] = z5i;
          ff[i__ + n5] = z6i;
          ff[i__ + n6] = z7i;
          zz[i__] = rcoeff7.t711 * z1i + rcoeff7.t712 * z2i +
                    rcoeff7.t713 * z3i + rcoeff7.t714 * z4i +
                    rcoeff7.t715 * z5i + rcoeff7.t716 * z6i +
                    rcoeff7.t717 * z7i;
          zz[i__ + n] = rcoeff7.t721 * z1i + rcoeff7.t722 * z2i +
                        rcoeff7.t723 * z3i + rcoeff7.t724 * z4i +
                        rcoeff7.t725 * z5i + rcoeff7.t726 * z6i +
                        rcoeff7.t727 * z7i;
          zz[i__ + n2] = rcoeff7.t731 * z1i + rcoeff7.t732 * z2i +
                         rcoeff7.t733 * z3i + rcoeff7.t734 * z4i +
                         rcoeff7.t735 * z5i + rcoeff7.t736 * z6i +
                         rcoeff7.t737 * z7i;
          zz[i__ + n3] = rcoeff7.t741 * z1i + rcoeff7.t742 * z2i +
                         rcoeff7.t743 * z3i + rcoeff7.t744 * z4i +
                         rcoeff7.t745 * z5i + rcoeff7.t746 * z6i +
                         rcoeff7.t747 * z7i;
          zz[i__ + n4] = rcoeff7.t751 * z1i + rcoeff7.t752 * z2i +
                         rcoeff7.t753 * z3i + rcoeff7.t754 * z4i +
                         rcoeff7.t755 * z5i + rcoeff7.t756 * z6i +
                         rcoeff7.t757 * z7i;
          zz[i__ + n5] = rcoeff7.t761 * z1i + rcoeff7.t762 * z2i +
                         rcoeff7.t763 * z3i + rcoeff7.t764 * z4i +
                         rcoeff7.t765 * z5i + rcoeff7.t766 * z6i +
                         rcoeff7.t767 * z7i;
          zz[i__ + n6] = rcoeff7.t771 * z1i + z2i + z4i + z6i;
        }
        if (faccon * dyno > fnewt) {
          goto L240;
        }
        /* --- ERROR ESTIMATION */
        estrav_(n, &fjac[fjac_offset], ldjac, mljac, mujac, &fmas[fmas_offset],
                ldmas, mlmas, mumas, h, dd, (S_fp)fcn, nfcn, &y0[1], &y[1],
                ijob, x, m1, m2, nm1, ns, &nns, &e1[e1_offset], lde1, &zz[1],
                &cont[1], &ff[1], &ip1[1], &iphes[1], &scal[1], &err, &first,
                &reject, &fac1, &rpar[1], &ipar[1]);
        /*       --- COMPUTE FINITE DIFFERENCES FOR DENSE OUTPUT */
        if (err < 1.) {
          i__1 = n;
          for (i__ = 1; i__ <= i__1; ++i__) {
            y[i__] += zz[i__ + (ns - 1) * n];
            cont[i__ + ns * n] = zz[i__] / weight_1.c__[1];
          }
          i__1 = ns - 1;
          for (k = 1; k <= i__1; ++k) {
            fact = 1. / (weight_1.c__[ns - k] - weight_1.c__[ns - k + 1]);
            i__2 = n;
            for (i__ = 1; i__ <= i__2; ++i__) {
              cont[i__ + k * n] =
                  (zz[i__ + (ns - k - 1) * n] - zz[i__ + (ns - k) * n]) * fact;
            }
          }
          i__1 = ns;
          for (j = 2; j <= i__1; ++j) {
            i__2 = j;
            for (k = ns; k >= i__2; --k) {
              fact = 1. / (weight_1.c__[ns - k] - weight_1.c__[ns - k + j]);
              i__3 = n;
              for (i__ = 1; i__ <= i__3; ++i__) {
                cont[i__ + k * n] =
                    (cont[i__ + k * n] - cont[i__ + (k - 1) * n]) * fact;
              }
            }
          }
        }
        /* *** *** *** *** *** *** *** */
        /* *** *** *** *** *** *** *** */
      } else {
        /* ASE       (NS.EQ.1) */
        /* *** *** *** *** *** *** *** */
        /* *** *** *** *** *** *** *** */
        if (first || *startn || change) {
          i__1 = ns;
          for (i__ = 1; i__ <= i__1; ++i__) {
            zz[i__] = 0.;
            ff[i__] = 0.;
          }
        } else {
          hquot = *h / hold;
          i__1 = n;
          for (i__ = 1; i__ <= i__1; ++i__) {
            z1i = hquot * cont[i__ + n];
            zz[i__] = z1i;
            ff[i__] = z1i;
          }
        }
        /* *** *** *** *** *** *** *** */
        /*  LOOP FOR THE SIMPLIFIED NEWTON ITERATION */
        /* *** *** *** *** *** *** *** */
        newt = 0;
        nit = nit1 - 3;
        expmi = 1. / expmns;
        /* Computing MAX */
        d__1 = *uround * 10 / rtol1;
        fnewt = max(d__1, .03);
        d__1 = max(faccon, *uround);
        faccon = pow_dd(&d__1, &c_b100);
        theta = abs(*thet);
      L440:
        if (newt >= nit) {
          goto L78;
        }
        /* ---     COMPUTE THE RIGHT-HAND SIDE */
        i__1 = n;
        for (i__ = 1; i__ <= i__1; ++i__) {
          cont[i__] = y[i__] + zz[i__];
        }
        (*fcn)(n, &xph, &cont[1], &zz[1], &rpar[1], &ipar[1]);
        ++(nfcn);
        /* ---     SOLVE THE LINEAR SYSTEMS */
        slvrar_(n, &fjac[fjac_offset], ldjac, mljac, mujac, &fmas[fmas_offset],
                ldmas, mlmas, mumas, m1, m2, nm1, &fac1, &e1[e1_offset], lde1,
                &zz[1], &ff[1], &ip1[1], &iphes[1], &ier, ijob);
        ++(nsol);
        ++newt;
        dyno = 0.;
        i__1 = n;
        for (i__ = 1; i__ <= i__1; ++i__) {
          denom = scal[i__];
          /* Computing 2nd power */
          d__1 = zz[i__] / denom;
          dyno += d__1 * d__1;
        }
        dyno = sqrt(dyno / nns);
        /* ---     BAD CONVERGENCE OR NUMBER OF ITERATIONS TO LARGE */
        if (newt > 1 && newt < nit) {
          thq = dyno / dynold;
          if (newt == 2) {
            theta = thq;
          } else {
            theta = sqrt(thq * thqold);
          }
          thqold = thq;
          if (theta < .99) {
            faccon = theta / (1. - theta);
            i__1 = nit - 1 - newt;
            dyth = faccon * dyno * pow_di(&theta, &i__1) / fnewt;
            if (dyth >= 1.) {
              /* Computing MAX */
              d__1 = 1e-4, d__2 = min(20., dyth);
              qnewt = max(d__1, d__2);
              d__1 = -1. / (nit + 4. - 1 - newt);
              hhfac = pow_dd(&qnewt, &d__1) * .8;
              *h = hhfac * *h;
              reject = true;
              last = false;
              if (hhfac <= .5) {
                unexn = true;
              }
              if (caljac) {
                goto L20;
              }
              goto L10;
            }
          } else {
            goto L78;
          }
        }
        dynold = max(dyno, *uround);
        i__1 = n;
        for (i__ = 1; i__ <= i__1; ++i__) {
          f1i = ff[i__] + zz[i__];
          ff[i__] = f1i;
          zz[i__] = f1i;
        }
        if (faccon * dyno > fnewt) {
          goto L440;
        }
        /* --- ERROR ESTIMATION */
        estrav_(n, &fjac[fjac_offset], ldjac, mljac, mujac, &fmas[fmas_offset],
                ldmas, mlmas, mumas, h, dd, (S_fp)fcn, nfcn, &y0[1], &y[1],
                ijob, x, m1, m2, nm1, ns, &nns, &e1[e1_offset], lde1, &zz[1],
                &cont[1], &ff[1], &ip1[1], &iphes[1], &scal[1], &err, &first,
                &reject, &fac1, &rpar[1], &ipar[1]);
        /*       --- COMPUTE FINITE DIFFERENCES FOR DENSE OUTPUT */
        if (err < 1.) {
          i__1 = n;
          for (i__ = 1; i__ <= i__1; ++i__) {
            y[i__] += zz[i__];
            cont[i__ + n] = zz[i__];
          }
        }
        /* *** *** *** *** *** *** *** */
        /* *** *** *** *** *** *** *** */
      }
    }
  }
  /* *** *** *** *** *** *** *** */
  /* *** *** *** *** *** *** *** */
  /* --- COMPUTATION OF HNEW */
  /* --- WE REQUIRE .2<=HNEW/H<=8. */
  /* Computing MIN */
  d__1 = *safe, d__2 = ((nit << 1) + 1) * *safe / (newt + (nit << 1));
  fac = min(d__1, d__2);
  /* Computing MAX */
  /* Computing MIN */
  d__3 = *facl, d__4 = pow_dd(&err, &expo) / fac;
  d__1 = *facr, d__2 = min(d__3, d__4);
  quot = max(d__1, d__2);
  hnew = *h / quot;
  /* *** *** *** *** *** *** *** */
  /*  IS THE ERROR SMALL ENOUGH ? */
  /* *** *** *** *** *** *** *** */
  if (err < 1.) {
    /* --- STEP IS ACCEPTED */
    first = false;
    ++(naccpt);
    if (*pred && !change) {
      /*       --- PREDICTIVE CONTROLLER OF GUSTAFSSON */
      if (naccpt > 1) {
        /* Computing 2nd power */
        d__2 = err;
        d__1 = d__2 * d__2 / erracc;
        facgus = hacc / *h * pow_dd(&d__1, &expo) / *safe;
        /* Computing MAX */
        d__1 = *facr, d__2 = min(*facl, facgus);
        facgus = max(d__1, d__2);
        quot = max(quot, facgus);
        hnew = *h / quot;
      }
      hacc = *h;
      erracc = max(.01, err);
    }
    xold = *x;
    hold = *h;
    *x = xph;
    /*       --- UPDATE SCALING */
    if (*itol == 0) {
      i__1 = n;
      for (i__ = 1; i__ <= i__1; ++i__) {
        scal[i__] = atol1 + rtol1 * (d__1 = y[i__], abs(d__1));
      }
    } else {
      i__1 = n;
      for (i__ = 1; i__ <= i__1; ++i__) {
        quott = atol[i__] / rtol[i__];
        rtol1 = pow_dd(&rtol[i__], &expmns) * .1;
        atol1 = rtol1 * quott;
        scal[i__] = atol1 + rtol1 * (d__1 = y[i__], abs(d__1));
      }
    }
    if (*iout != 0) {
      nrsol = naccpt + 1;
      weight_1.xsol = *x;
      xosol = xold;
      i__1 = n;
      for (i__ = 1; i__ <= i__1; ++i__) {
        cont[i__] = y[i__];
      }
      nsolu = n;
      weight_1.hsol = hold;
      (*solout)(&nrsol, &xosol, &weight_1.xsol, &y[1], &cont[1], &lrc, &nsolu,
                &rpar[1], &ipar[1], &irtrn);
      if (irtrn < 0) {
        goto L179;
      }
    }
    caljac = false;
    if (last) {
      *h = hopt;
      *idid = 1;
      return 0;
    }
    (*fcn)(n, x, &y[1], &y0[1], &rpar[1], &ipar[1]);
    ++(nfcn);
    /* Computing MIN */
    d__1 = abs(hnew);
    hnew = posneg * min(d__1, hmaxn);
    hopt = hnew;
    hopt = min(*h, hnew);
    if (reject) {
      /* Computing MIN */
      d__1 = abs(hnew), d__2 = abs(*h);
      hnew = posneg * min(d__1, d__2);
    }
    reject = false;
    if ((*x + hnew / *quot1 - *xend) * posneg >= 0.) {
      *h = *xend - *x;
      last = true;
    } else {
      qt = hnew / *h;
      hhfac = *h;
      if (theta <= *thet && qt >= *quot1 && qt <= *quot2) {
        ikeep = 1;
        goto L30;
      }
      *h = hnew;
    }
    hhfac = *h;
    if (theta <= *thet) {
      goto L20;
    }
    goto L10;
  } else {
    /* --- STEP IS REJECTED */
    reject = true;
    last = false;
    if (first) {
      *h *= .1;
      hhfac = .1;
    } else {
      hhfac = hnew / *h;
      *h = hnew;
    }
    if (naccpt >= 1) {
      ++(nrejct);
    }
    if (caljac) {
      goto L20;
    }
    goto L10;
  }
/* --- UNEXPECTED STEP-REJECTION */
L78:
  unexp = true;
  if (ier != 0) {
    ++nsing;
    if (nsing >= 5) {
      goto L176;
    }
  }
  *h *= .5;
  hhfac = .5;
  reject = true;
  last = false;
  if (caljac) {
    goto L20;
  }
  goto L10;
/* --- FAIL EXIT */
L176:
  s_wsfe(&io___195);
  do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(double));
  e_wsfe();
  s_wsle(&io___196);
  do_lio(&c__9, &c__1, " MATRIX IS REPEATEDLY SINGULAR, IER=", (ftnlen)36);
  do_lio(&c__3, &c__1, (char *)&ier, (ftnlen)sizeof(int));
  e_wsle();
  *idid = -4;
  return 0;
L177:
  s_wsfe(&io___197);
  do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(double));
  e_wsfe();
  s_wsle(&io___198);
  do_lio(&c__9, &c__1, " STEP SIZE T0O SMALL, H=", (ftnlen)24);
  do_lio(&c__5, &c__1, (char *)&(*h), (ftnlen)sizeof(double));
  e_wsle();
  *idid = -3;
  return 0;
L178:
  s_wsfe(&io___199);
  do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(double));
  e_wsfe();
  s_wsle(&io___200);
  do_lio(&c__9, &c__1, " MORE THAN NMAX =", (ftnlen)17);
  do_lio(&c__3, &c__1, (char *)&(nmax), (ftnlen)sizeof(int));
  do_lio(&c__9, &c__1, "STEPS ARE NEEDED", (ftnlen)16);
  e_wsle();
  *idid = -2;
  return 0;
/* --- EXIT CAUSED BY SOLOUT */
L179:
  s_wsfe(&io___201);
  do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(double));
  e_wsfe();
  *idid = 2;
  return 0;
}

int radau(int n, U_fp fcn, double *x, double *y, double *xend, double *h__,
          double *rtol, double *atol, int *itol, U_fp jac, int *ijac,
          int *mljac, int *mujac, U_fp mas, int *imas, int *mlmas, int *mumas,
          U_fp solout, int *iout, double *work, int *lwork, int *iwork,
          int *liwork, double *rpar, int *ipar, int *idid) {
  /* System generated locals */
  int i__1;

  /* Builtin functions */
  int s_wsle(cilist *), do_lio(int *, int *, char *, ftnlen), e_wsle();

  /* Local variables */
  int i__, m1, m2, ns, nm1, iee, nit, nns, iee1, lde1, iey0;
  double facl;
  int ndec, njac;
  double facr;
  int ieff;
  double safe, hhod;
  int ijob, nfcn, nmee;
  bool pred;
  double hmax;
  int nmax;
  double thet, vitd, hhou;
  int nsol, iezz;
  double vitu;
  int nsus, ieip1, ieip2, nind1, nind2, nind3, nm1ns;
  double quot1, quot2;
  int iejac, ldjac;
  bool jband;
  int iecon, iemas, ldmas, ieiph;
  bool arret;
  int nsmin, nsmax, nstep, ldmas2, iescal, naccpt;
  int nrejct;
  bool implct;
  int istore;
  bool startn;
  double uround;

  /* Fortran I/O blocks */
  cilist io___13 = {0, 6, 0, 0, 0};
  cilist io___15 = {0, 6, 0, 0, 0};
  cilist io___17 = {0, 6, 0, 0, 0};
  cilist io___22 = {0, 6, 0, 0, 0};
  cilist io___27 = {0, 6, 0, 0, 0};
  cilist io___29 = {0, 6, 0, 0, 0};
  cilist io___30 = {0, 6, 0, 0, 0};
  cilist io___32 = {0, 6, 0, 0, 0};
  cilist io___34 = {0, 6, 0, 0, 0};
  cilist io___36 = {0, 6, 0, 0, 0};
  cilist io___39 = {0, 6, 0, 0, 0};
  cilist io___43 = {0, 6, 0, 0, 0};
  cilist io___54 = {0, 6, 0, 0, 0};
  cilist io___56 = {0, 6, 0, 0, 0};
  cilist io___70 = {0, 6, 0, 0, 0};
  cilist io___74 = {0, 6, 0, 0, 0};

  /* ---------------------------------------------------------- */
  /*     NUMERICAL SOLUTION OF A STIFF (OR DIFFERENTIAL ALGEBRAIC) */
  /*     SYSTEM OF FIRST 0RDER ORDINARY DIFFERENTIAL EQUATIONS */
  /*                     M*Y'=F(X,Y). */
  /*     THE SYSTEM CAN BE (LINEARLY) IMPLICIT (MASS-MATRIX M .NE. I) */
  /*     OR EXPLICIT (M=I). */
  /*     THE CODE IS BASED ON IMPLICIT RUNGE-KUTTA METHODS (RADAU IIA) */
  /*     WITH VARIABLE ORDER (5, 9, 13), WITH STEP SIZE CONTROL */
  /*     AND CONTINUOUS OUTPUT. */

  /*     AUTHORS: E. HAIRER AND G. WANNER */
  /*              UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES */
  /*              CH-1211 GENEVE 24, SWITZERLAND */
  /*              E-MAIL:  Ernst.Hairer@math.unige.ch */
  /*                       Gerhard.Wanner@math.unige.ch */

  /*     FOR A DESCRIPTION OF THE RELATED CODE RADAU5 SEE THE BOOK: */
  /*         E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL */
  /*         EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS. */
  /*         SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS 14, */
  /*         SPRINGER-VERLAG 1991, SECOND EDITION 1996. */

  /*     PRELIMINARY VERSION OF APRIL 23, 1998 */
  /*     (latest small correction: January 18, 2002) */

  /*     INPUT PARAMETERS */
  /*     ---------------- */
  /*     N           DIMENSION OF THE SYSTEM */

  /*     FCN         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE */
  /*                 VALUE OF F(X,Y): */
  /*                    SUBROUTINE FCN(N,X,Y,F,RPAR,IPAR) */
  /*                    DOUBLE PRECISION X,Y(N),F(N) */
  /*                    F(1)=...   ETC. */
  /*                 RPAR, IPAR (SEE BELOW) */

  /*     X           INITIAL X-VALUE */

  /*     Y(N)        INITIAL VALUES FOR Y */

  /*     XEND        FINAL X-VALUE (XEND-X MAY BE POSITIVE OR NEGATIVE) */

  /*     H           INITIAL STEP SIZE GUESS; */
  /*                 FOR STIFF EQUATIONS WITH INITIAL TRANSIENT, */
  /*                 H=1.D0/(NORM OF F'), USUALLY 1.D-3 OR 1.D-5, IS GOOD. */
  /*                 THIS CHOICE IS NOT VERY IMPORTANT, THE STEP SIZE IS */
  /*                 QUICKLY ADAPTED. (IF H=0.D0, THE CODE PUTS H=1.D-6). */

  /*     RTOL,ATOL   RELATIVE AND ABSOLUTE ERROR TOLERANCES. THEY */
  /*                 CAN BE BOTH SCALARS OR ELSE BOTH VECTORS OF LENGTH N. */

  /*     ITOL        SWITCH FOR RTOL AND ATOL: */
  /*                   ITOL=0: BOTH RTOL AND ATOL ARE SCALARS. */
  /*                     THE CODE KEEPS, ROUGHLY, THE LOCAL ERROR OF */
  /*                     Y(I) BELOW RTOL*ABS(Y(I))+ATOL */
  /*                   ITOL=1: BOTH RTOL AND ATOL ARE VECTORS. */
  /*                     THE CODE KEEPS THE LOCAL ERROR OF Y(I) BELOW */
  /*                     RTOL(I)*ABS(Y(I))+ATOL(I). */

  /*     JAC         NAME (EXTERNAL) OF THE SUBROUTINE WHICH COMPUTES */
  /*                 THE PARTIAL DERIVATIVES OF F(X,Y) WITH RESPECT TO Y */
  /*                 (THIS ROUTINE IS ONLY CALLED IF IJAC=1; SUPPLY */
  /*                 A DUMMY SUBROUTINE IN THE CASE IJAC=0). */
  /*                 FOR IJAC=1, THIS SUBROUTINE MUST HAVE THE FORM */
  /*                    SUBROUTINE JAC(N,X,Y,DFY,LDFY,RPAR,IPAR) */
  /*                    DOUBLE PRECISION X,Y(N),DFY(LDFY,N) */
  /*                    DFY(1,1)= ... */
  /*                 LDFY, THE COLUMN-LENGTH OF THE ARRAY, IS */
  /*                 FURNISHED BY THE CALLING PROGRAM. */
  /*                 IF (MLJAC.EQ.N) THE JACOBIAN IS SUPPOSED TO */
  /*                    BE FULL AND THE PARTIAL DERIVATIVES ARE */
  /*                    STORED IN DFY AS */
  /*                       DFY(I,J) = PARTIAL F(I) / PARTIAL Y(J) */
  /*                 ELSE, THE JACOBIAN IS TAKEN AS BANDED AND */
  /*                    THE PARTIAL DERIVATIVES ARE STORED */
  /*                    DIAGONAL-WISE AS */
  /*                       DFY(I-J+MUJAC+1,J) = PARTIAL F(I) / PARTIAL Y(J). */

  /*     IJAC        SWITCH FOR THE COMPUTATION OF THE JACOBIAN: */
  /*                    IJAC=0: JACOBIAN IS COMPUTED INTERNALLY BY FINITE */
  /*                       DIFFERENCES, SUBROUTINE "JAC" IS NEVER CALLED. */
  /*                    IJAC=1: JACOBIAN IS SUPPLIED BY SUBROUTINE JAC. */

  /*     MLJAC       SWITCH FOR THE BANDED STRUCTURE OF THE JACOBIAN: */
  /*                    MLJAC=N: JACOBIAN IS A FULL MATRIX. THE LINEAR */
  /*                       ALGEBRA IS DONE BY FULL-MATRIX GAUSS-ELIMINATION. */
  /*                    0<=MLJAC<N: MLJAC IS THE LOWER BANDWITH OF JACOBIAN */
  /*                       MATRIX (>= NUMBER OF NON-ZERO DIAGONALS BELOW */
  /*                       THE MAIN DIAGONAL). */

  /*     MUJAC       UPPER BANDWITH OF JACOBIAN  MATRIX (>= NUMBER OF NON- */
  /*                 ZERO DIAGONALS ABOVE THE MAIN DIAGONAL). */
  /*                 NEED NOT BE DEFINED IF MLJAC=N. */

  /*     ----   MAS,IMAS,MLMAS, AND MUMAS HAVE ANALOG MEANINGS      ----- */
  /*     ----   FOR THE "MASS MATRIX" (THE MATRIX "M" OF SECTION IV.8): - */

  /*     MAS         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE MASS- */
  /*                 MATRIX M. */
  /*                 IF IMAS=0, THIS MATRIX IS ASSUMED TO BE THE IDENTITY */
  /*                 MATRIX AND NEEDS NOT TO BE DEFINED; */
  /*                 SUPPLY A DUMMY SUBROUTINE IN THIS CASE. */
  /*                 IF IMAS=1, THE SUBROUTINE MAS IS OF THE FORM */
  /*                    SUBROUTINE MAS(N,AM,LMAS,RPAR,IPAR) */
  /*                    DOUBLE PRECISION AM(LMAS,N) */
  /*                    AM(1,1)= .... */
  /*                    IF (MLMAS.EQ.N) THE MASS-MATRIX IS STORED */
  /*                    AS FULL MATRIX LIKE */
  /*                         AM(I,J) = M(I,J) */
  /*                    ELSE, THE MATRIX IS TAKEN AS BANDED AND STORED */
  /*                    DIAGONAL-WISE AS */
  /*                         AM(I-J+MUMAS+1,J) = M(I,J). */

  /*     IMAS       GIVES INFORMATION ON THE MASS-MATRIX: */
  /*                    IMAS=0: M IS SUPPOSED TO BE THE IDENTITY */
  /*                       MATRIX, MAS IS NEVER CALLED. */
  /*                    IMAS=1: MASS-MATRIX  IS SUPPLIED. */

  /*     MLMAS       SWITCH FOR THE BANDED STRUCTURE OF THE MASS-MATRIX: */
  /*                    MLMAS=N: THE FULL MATRIX CASE. THE LINEAR */
  /*                       ALGEBRA IS DONE BY FULL-MATRIX GAUSS-ELIMINATION. */
  /*                    0<=MLMAS<N: MLMAS IS THE LOWER BANDWITH OF THE */
  /*                       MATRIX (>= NUMBER OF NON-ZERO DIAGONALS BELOW */
  /*                       THE MAIN DIAGONAL). */
  /*                 MLMAS IS SUPPOSED TO BE .LE. MLJAC. */

  /*     MUMAS       UPPER BANDWITH OF MASS-MATRIX (>= NUMBER OF NON- */
  /*                 ZERO DIAGONALS ABOVE THE MAIN DIAGONAL). */
  /*                 NEED NOT BE DEFINED IF MLMAS=N. */
  /*                 MUMAS IS SUPPOSED TO BE .LE. MUJAC. */

  /*     SOLOUT      NAME (EXTERNAL) OF SUBROUTINE PROVIDING THE */
  /*                 NUMERICAL SOLUTION DURING INTEGRATION. */
  /*                 IF IOUT=1, IT IS CALLED AFTER EVERY SUCCESSFUL STEP. */
  /*                 SUPPLY A DUMMY SUBROUTINE IF IOUT=0. */
  /*                 IT MUST HAVE THE FORM */
  /*                    SUBROUTINE SOLOUT (NR,XOLD,X,Y,CONT,LRC,N, */
  /*                                       RPAR,IPAR,IRTRN) */
  /*                    DOUBLE PRECISION X,Y(N),CONT(LRC) */
  /*                    .... */
  /*                 SOLOUT FURNISHES THE SOLUTION "Y" AT THE NR-TH */
  /*                    GRID-POINT "X" (THEREBY THE INITIAL VALUE IS */
  /*                    THE FIRST GRID-POINT). */
  /*                 "XOLD" IS THE PRECEEDING GRID-POINT. */
  /*                 "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN */
  /*                    IS SET <0, RADAU RETURNS TO THE CALLING PROGRAM. */

  /*          -----  CONTINUOUS OUTPUT: ----- */
  /*                 DURING CALLS TO "SOLOUT", A CONTINUOUS SOLUTION */
  /*                 FOR THE INTERVAL [XOLD,X] IS AVAILABLE THROUGH */
  /*                 THE FUNCTION */
  /*                        >>>   CONTRA(I,S,CONT,LRC)   <<< */
  /*                 WHICH PROVIDES AN APPROXIMATION TO THE I-TH */
  /*                 COMPONENT OF THE SOLUTION AT THE POINT S. THE VALUE */
  /*                 S SHOULD LIE IN THE INTERVAL [XOLD,X]. */

  /*     IOUT        SWITCH FOR CALLING THE SUBROUTINE SOLOUT: */
  /*                    IOUT=0: SUBROUTINE IS NEVER CALLED */
  /*                    IOUT=1: SUBROUTINE IS AVAILABLE FOR OUTPUT. */

  /*     WORK        ARRAY OF WORKING SPACE OF LENGTH "LWORK". */
  /*                 WORK(1), WORK(2),.., WORK(20) SERVE AS PARAMETERS */
  /*                 FOR THE CODE. FOR STANDARD USE OF THE CODE */
  /*                 WORK(1),..,WORK(20) MUST BE SET TO ZERO BEFORE THE */
  /*                 FIRST CALL. SEE BELOW FOR A MORE SOPHISTICATED USE. */
  /*                 WORK(8),..,WORK(LWORK) SERVE AS WORKING SPACE */
  /*                 FOR ALL VECTORS AND MATRICES. */
  /*                 "LWORK" MUST BE AT LEAST */
  /*                          N*(LJAC+LMAS+NSMAX*LE+3*NSMAX+3)+20 */
  /*                 WHERE */
  /*                    NSMAX=IWORK(12) (SEE BELOW) */
  /*                 AND */
  /*                    LJAC=N              IF MLJAC=N (FULL JACOBIAN) */
  /*                    LJAC=MLJAC+MUJAC+1  IF MLJAC<N (BANDED JAC.) */
  /*                 AND */
  /*                    LMAS=0              IF IMAS=0 */
  /*                    LMAS=N              IF IMAS=1 AND MLMAS=N (FULL) */
  /*                    LMAS=MLMAS+MUMAS+1  IF MLMAS<N (BANDED MASS-M.) */
  /*                 AND */
  /*                    LE=N               IF MLJAC=N (FULL JACOBIAN) */
  /*                    LE=2*MLJAC+MUJAC+1 IF MLJAC<N (BANDED JAC.) */

  /*                 IN THE USUAL CASE WHERE THE JACOBIAN IS FULL AND THE */
  /*                 MASS-MATRIX IS THE INDENTITY (IMAS=0), THE MINIMUM */
  /*                 STORAGE REQUIREMENT IS */
  /*                      LWORK = (NSMAX+1)*N*N+(3*NSMAX+3)*N+20. */
  /*                 IF IWORK(9)=M1>0 THEN "LWORK" MUST BE AT LEAST */
  /*                      N*(LJAC+3*NSMAX+3)+(N-M1)*(LMAS+NSMAX*LE)+20 */
  /*                 WHERE IN THE DEFINITIONS OF LJAC, LMAS AND LE THE */
  /*                 NUMBER N CAN BE REPLACED BY N-M1. */

  /*     LWORK       DECLARED LENGTH OF ARRAY "WORK". */

  /*     IWORK       INTEGER WORKING SPACE OF LENGTH "LIWORK". */
  /*                 IWORK(1),IWORK(2),...,IWORK(20) SERVE AS PARAMETERS */
  /*                 FOR THE CODE. FOR STANDARD USE, SET IWORK(1),.., */
  /*                 IWORK(20) TO ZERO BEFORE CALLING. */
  /*                 IWORK(21),...,IWORK(LIWORK) SERVE AS WORKING AREA. */
  /*                 "LIWORK" MUST BE AT LEAST */
  /*                             (2+(NSMAX-1)/2)*N+20. */

  /*     LIWORK      DECLARED LENGTH OF ARRAY "IWORK". */

  /*     RPAR, IPAR  REAL AND INTEGER PARAMETERS (OR PARAMETER ARRAYS) WHICH */
  /*                 CAN BE USED FOR COMMUNICATION BETWEEN YOUR CALLING */
  /*                 PROGRAM AND THE FCN, JAC, MAS, SOLOUT SUBROUTINES. */

  /* ---------------------------------------------------------------------- */

  /*     SOPHISTICATED SETTING OF PARAMETERS */
  /*     ----------------------------------- */
  /*              SEVERAL PARAMETERS OF THE CODE ARE TUNED TO MAKE IT WORK */
  /*              WELL. THEY MAY BE DEFINED BY SETTING WORK(1),... */
  /*              AS WELL AS IWORK(1),... DIFFERENT FROM ZERO. */
  /*              FOR ZERO INPUT, THE CODE CHOOSES DEFAULT VALUES: */

  /*    IWORK(1)  IF IWORK(1).NE.0, THE CODE TRANSFORMS THE JACOBIAN */
  /*              MATRIX TO HESSENBERG FORM. THIS IS PARTICULARLY */
  /*              ADVANTAGEOUS FOR LARGE SYSTEMS WITH FULL JACOBIAN. */
  /*              IT DOES NOT WORK FOR BANDED JACOBIAN (MLJAC<N) */
  /*              AND NOT FOR IMPLICIT SYSTEMS (IMAS=1). */

  /*    IWORK(2)  THIS IS THE MAXIMAL NUMBER OF ALLOWED STEPS. */
  /*              THE DEFAULT VALUE (FOR IWORK(2)=0) IS 100000. */

  /*    IWORK(3)  THE MAXIMUM NUMBER OF NEWTON ITERATIONS FOR THE */
  /*              SOLUTION OF THE IMPLICIT SYSTEM IN EACH STEP */
  /*              IWORK(3)+(NS-3)*2.5. DEFAULT VALUE (FOR IWORK(3)=0) IS 7. */
  /*              NS IS THE NUMBER OF STAGES (SEE IWORK(11)). */

  /*    IWORK(4)  IF IWORK(4).EQ.0 THE EXTRAPOLATED COLLOCATION SOLUTION */
  /*              IS TAKEN AS STARTING VALUE FOR NEWTON'S METHOD. */
  /*              IF IWORK(4).NE.0 ZERO STARTING VALUES ARE USED. */
  /*              THE LATTER IS RECOMMENDED IF NEWTON'S METHOD HAS */
  /*              DIFFICULTIES WITH CONVERGENCE (THIS IS THE CASE WHEN */
  /*              NSTEP IS LARGER THAN NACCPT + NREJCT; SEE OUTPUT PARAM.). */
  /*              DEFAULT IS IWORK(4)=0. */

  /*       THE FOLLOWING 3 PARAMETERS ARE IMPORTANT FOR */
  /*       DIFFERENTIAL-ALGEBRAIC SYSTEMS OF INDEX > 1. */
  /*       THE FUNCTION-SUBROUTINE SHOULD BE WRITTEN SUCH THAT */
  /*       THE INDEX 1,2,3 VARIABLES APPEAR IN THIS ORDER. */
  /*       IN ESTIMATING THE ERROR THE INDEX 2 VARIABLES ARE */
  /*       MULTIPLIED BY H, THE INDEX 3 VARIABLES BY H**2. */

  /*    IWORK(5)  DIMENSION OF THE INDEX 1 VARIABLES (MUST BE > 0). FOR */
  /*              ODE'S THIS EQUALS THE DIMENSION OF THE SYSTEM. */
  /*              DEFAULT IWORK(5)=N. */

  /*    IWORK(6)  DIMENSION OF THE INDEX 2 VARIABLES. DEFAULT IWORK(6)=0. */

  /*    IWORK(7)  DIMENSION OF THE INDEX 3 VARIABLES. DEFAULT IWORK(7)=0. */

  /*    IWORK(8)  SWITCH FOR STEP SIZE STRATEGY */
  /*              IF IWORK(8).EQ.1  MOD. PREDICTIVE CONTROLLER (GUSTAFSSON) */
  /*              IF IWORK(8).EQ.2  CLASSICAL STEP SIZE CONTROL */
  /*              THE DEFAULT VALUE (FOR IWORK(8)=0) IS IWORK(8)=1. */
  /*              THE CHOICE IWORK(8).EQ.1 SEEMS TO PRODUCE SAFER RESULTS; */
  /*              FOR SIMPLE PROBLEMS, THE CHOICE IWORK(8).EQ.2 PRODUCES */
  /*              OFTEN SLIGHTLY FASTER RUNS */

  /*       IF THE DIFFERENTIAL SYSTEM HAS THE SPECIAL STRUCTURE THAT */
  /*            Y(I)' = Y(I+M2)   FOR  I=1,...,M1, */
  /*       WITH M1 A MULTIPLE OF M2, A SUBSTANTIAL GAIN IN COMPUTERTIME */
  /*       CAN BE ACHIEVED BY SETTING THE PARAMETERS IWORK(9) AND IWORK(10). */
  /*       E.G., FOR SECOND ORDER SYSTEMS P'=V, V'=G(P,V), WHERE P AND V ARE */
  /*       VECTORS OF DIMENSION N/2, ONE HAS TO PUT M1=M2=N/2. */
  /*       FOR M1>0 SOME OF THE INPUT PARAMETERS HAVE DIFFERENT MEANINGS: */
  /*       - JAC: ONLY THE ELEMENTS OF THE NON-TRIVIAL PART OF THE */
  /*              JACOBIAN HAVE TO BE STORED */
  /*              IF (MLJAC.EQ.N-M1) THE JACOBIAN IS SUPPOSED TO BE FULL */
  /*                 DFY(I,J) = PARTIAL F(I+M1) / PARTIAL Y(J) */
  /*                FOR I=1,N-M1 AND J=1,N. */
  /*              ELSE, THE JACOBIAN IS BANDED ( M1 = M2 * MM ) */
  /*                 DFY(I-J+MUJAC+1,J+K*M2) = PARTIAL F(I+M1) / PARTIAL
   * Y(J+K*M2) */
  /*                FOR I=1,MLJAC+MUJAC+1 AND J=1,M2 AND K=0,MM. */
  /*       - MLJAC: MLJAC=N-M1: IF THE NON-TRIVIAL PART OF THE JACOBIAN IS FULL
   */
  /*                0<=MLJAC<N-M1: IF THE (MM+1) SUBMATRICES (FOR K=0,MM) */
  /*                     PARTIAL F(I+M1) / PARTIAL Y(J+K*M2),  I,J=1,M2 */
  /*                    ARE BANDED, MLJAC IS THE MAXIMAL LOWER BANDWIDTH */
  /*                    OF THESE MM+1 SUBMATRICES */
  /*       - MUJAC: MAXIMAL UPPER BANDWIDTH OF THESE MM+1 SUBMATRICES */
  /*                NEED NOT BE DEFINED IF MLJAC=N-M1 */
  /*       - MAS: IF IMAS=0 THIS MATRIX IS ASSUMED TO BE THE IDENTITY AND */
  /*              NEED NOT BE DEFINED. SUPPLY A DUMMY SUBROUTINE IN THIS CASE.
   */
  /*              IT IS ASSUMED THAT ONLY THE ELEMENTS OF RIGHT LOWER BLOCK OF
   */
  /*              DIMENSION N-M1 DIFFER FROM THAT OF THE IDENTITY MATRIX. */
  /*              IF (MLMAS.EQ.N-M1) THIS SUBMATRIX IS SUPPOSED TO BE FULL */
  /*                 AM(I,J) = M(I+M1,J+M1)     FOR I=1,N-M1 AND J=1,N-M1. */
  /*              ELSE, THE MASS MATRIX IS BANDED */
  /*                 AM(I-J+MUMAS+1,J) = M(I+M1,J+M1) */
  /*       - MLMAS: MLMAS=N-M1: IF THE NON-TRIVIAL PART OF M IS FULL */
  /*                0<=MLMAS<N-M1: LOWER BANDWIDTH OF THE MASS MATRIX */
  /*       - MUMAS: UPPER BANDWIDTH OF THE MASS MATRIX */
  /*                NEED NOT BE DEFINED IF MLMAS=N-M1 */

  /*    IWORK(9)  THE VALUE OF M1.  DEFAULT M1=0. */

  /*    IWORK(10) THE VALUE OF M2.  DEFAULT M2=M1. */

  /*    IWORK(11) NSMIN, MINIMAL NUMBER OF STAGES NS (ORDER 2*NS-1) */
  /*              POSSIBLE VALUES ARE 1,3,5,7. DEFAULT NS=3. */

  /*    IWORK(12) NSMAX, MAXIMAL NUMBER OF STAGES NS. */
  /*              POSSIBLE VALUES ARE 1,3,5,7. DEFAULT NS=7. */

  /*    IWORK(13) VALUE OF NS FOR THE FIRST STEP (DEFAULT VALUE: NSMIN) */

  /* ---------- */

  /*    WORK(1)   UROUND, THE ROUNDING UNIT, DEFAULT 1.D-16. */

  /*    WORK(2)   THE SAFETY FACTOR IN STEP SIZE PREDICTION, */
  /*              DEFAULT 0.9D0. */

  /*    WORK(3)   DECIDES WHETHER THE JACOBIAN SHOULD BE RECOMPUTED; */
  /*              INCREASE WORK(3), TO 0.1 SAY, WHEN JACOBIAN EVALUATIONS */
  /*              ARE COSTLY. FOR SMALL SYSTEMS WORK(3) SHOULD BE SMALLER */
  /*              (0.001D0, SAY). NEGATIV WORK(3) FORCES THE CODE TO */
  /*              COMPUTE THE JACOBIAN AFTER EVERY ACCEPTED STEP. */
  /*              DEFAULT 0.001D0. */

  /*    WORK(5) AND WORK(6) : IF WORK(5) < HNEW/HOLD < WORK(6), THEN THE */
  /*              STEP SIZE IS NOT CHANGED. THIS SAVES, TOGETHER WITH A */
  /*              LARGE WORK(3), LU-DECOMPOSITIONS AND COMPUTING TIME FOR */
  /*              LARGE SYSTEMS. FOR SMALL SYSTEMS ONE MAY HAVE */
  /*              WORK(5)=1.D0, WORK(6)=1.2D0, FOR LARGE FULL SYSTEMS */
  /*              WORK(5)=0.99D0, WORK(6)=2.D0 MIGHT BE GOOD. */
  /*              DEFAULTS WORK(5)=1.D0, WORK(6)=1.2D0 . */

  /*    WORK(7)   MAXIMAL STEP SIZE, DEFAULT XEND-X. */

  /*    WORK(8), WORK(9)   PARAMETERS FOR STEP SIZE SELECTION */
  /*              THE NEW STEP SIZE IS CHOSEN SUBJECT TO THE RESTRICTION */
  /*                 WORK(8) <= HNEW/HOLD <= WORK(9) */
  /*              DEFAULT VALUES: WORK(8)=0.2D0, WORK(9)=8.D0 */

  /*    WORK(10)  ORDER IS INCREASED IF THE CONTRACTIVITY FACTOR IS */
  /*              SMALL THAN WORK(10), DEFAULT VALUE IS 0.002 */

  /*    WORK(11)  ORDER IS DECREASED IF THE CONTRACTIVITY FACTOR IS */
  /*              LARGER THAN WORK(11), DEFAULT VALUE IS 0.8 */

  /*    WORK(12), WORK(13)  ORDER IS DECREASED ONLY IF THE STEPSIZE */
  /*              RATIO SATISFIES  WORK(13).LE.HNEW/H.LE.WORK(12), */
  /*              DEFAULT VALUES ARE 1.2 AND 0.8 */

  /* ----------------------------------------------------------------------- */

  /*     OUTPUT PARAMETERS */
  /*     ----------------- */
  /*     X           X-VALUE FOR WHICH THE SOLUTION HAS BEEN COMPUTED */
  /*                 (AFTER SUCCESSFUL RETURN X=XEND). */

  /*     Y(N)        NUMERICAL SOLUTION AT X */

  /*     H           PREDICTED STEP SIZE OF THE LAST ACCEPTED STEP */

  /*     IDID        REPORTS ON SUCCESSFULNESS UPON RETURN: */
  /*                   IDID= 1  COMPUTATION SUCCESSFUL, */
  /*                   IDID= 2  COMPUT. SUCCESSFUL (INTERRUPTED BY SOLOUT) */
  /*                   IDID=-1  INPUT IS NOT CONSISTENT, */
  /*                   IDID=-2  LARGER NMAX IS NEEDED, */
  /*                   IDID=-3  STEP SIZE BECOMES TOO SMALL, */
  /*                   IDID=-4  MATRIX IS REPEATEDLY SINGULAR. */

  /*   IWORK(14)  NFCN    NUMBER OF FUNCTION EVALUATIONS (THOSE FOR NUMERICAL */
  /*                      EVALUATION OF THE JACOBIAN ARE NOT COUNTED) */
  /*   IWORK(15)  NJAC    NUMBER OF JACOBIAN EVALUATIONS (EITHER ANALYTICALLY */
  /*                      OR NUMERICALLY) */
  /*   IWORK(16)  NSTEP   NUMBER OF COMPUTED STEPS */
  /*   IWORK(17)  NACCPT  NUMBER OF ACCEPTED STEPS */
  /*   IWORK(18)  NREJCT  NUMBER OF REJECTED STEPS (DUE TO ERROR TEST), */
  /*                      (STEP REJECTIONS IN THE FIRST STEP ARE NOT COUNTED) */
  /*   IWORK(19)  NDEC    NUMBER OF LU-DECOMPOSITIONS OF THE MATRICES */
  /*   IWORK(20)  NSOL    NUMBER OF FORWARD-BACKWARD SUBSTITUTIONS, OF ALL */
  /*                      SYSTEMS THAT HAVE TO BE SOLVED FOR ONE SIMPLIFIED */
  /*                      NEWTON ITERATION; THE NSTEP (REAL) FORWARD-BACKWARD */
  /*                      SUBSTITUTIONS, NEEDED FOR STEP SIZE SELECTION, */
  /*                      ARE NOT COUNTED. */
  /* ----------------------------------------------------------------------- */
  /* *** *** *** *** *** *** *** *** *** *** *** *** *** */
  /*          DECLARATIONS */
  /* *** *** *** *** *** *** *** *** *** *** *** *** *** */
  /* *** *** *** *** *** *** *** */
  /*        SETTING THE PARAMETERS */
  /* *** *** *** *** *** *** *** */
  /* Parameter adjustments */
  --y;
  --rtol;
  --atol;
  --work;
  --iwork;
  --rpar;
  --ipar;

  /* Function Body */
  nfcn = 0;
  njac = 0;
  nstep = 0;
  naccpt = 0;
  nrejct = 0;
  ndec = 0;
  nsol = 0;
  arret = false;
  /* -------- NUMBER MAXIMAL AND MINIMAL OF STAGES  NS */
  if (iwork[11] == 0) {
    nsmin = 3;
  } else {
    nsmin = max(1, iwork[11]);
    if (iwork[11] >= 2) {
      nsmin = max(3, iwork[11]);
    }
    if (iwork[11] >= 4) {
      nsmin = max(5, iwork[11]);
    }
    if (iwork[11] >= 6) {
      nsmin = 7;
    }
  }
  if (iwork[12] == 0) {
    nsmax = 7;
  } else {
    nsmax = min(7, iwork[12]);
    if (iwork[12] <= 6) {
      nsmax = min(5, iwork[12]);
    }
    if (iwork[12] <= 4) {
      nsmax = min(3, iwork[12]);
    }
    if (iwork[12] <= 2) {
      nsmax = 1;
    }
  }
  ns = nsmax;
  if (iwork[13] == 0) {
    nsus = nsmin;
  } else {
    nsus = iwork[13];
    if (nsus <= 0 || ns >= 8 || ns == 2 || ns == 4 || ns == 6) {
      s_wsle(&io___13);
      do_lio(&c__9, &c__1, " WRONG INPUT IWORK(13)=", (ftnlen)23);
      do_lio(&c__3, &c__1, (char *)&iwork[13], (ftnlen)sizeof(int));
      e_wsle();
      arret = true;
    }
  }
  /* -------- NMAX , THE MAXIMAL NUMBER OF STEPS ----- */
  if (iwork[2] == 0) {
    nmax = 100000;
  } else {
    nmax = iwork[2];
    if (nmax <= 0) {
      s_wsle(&io___15);
      do_lio(&c__9, &c__1, " WRONG INPUT IWORK(2)=", (ftnlen)22);
      do_lio(&c__3, &c__1, (char *)&iwork[2], (ftnlen)sizeof(int));
      e_wsle();
      arret = true;
    }
  }
  /* -------- NIT    MAXIMAL NUMBER OF NEWTON ITERATIONS */
  if (iwork[3] == 0) {
    nit = 7;
  } else {
    nit = iwork[3];
    if (nit <= 0 || nit > 50) {
      s_wsle(&io___17);
      do_lio(&c__9, &c__1, " CURIOUS INPUT IWORK(3)=", (ftnlen)24);
      do_lio(&c__3, &c__1, (char *)&iwork[3], (ftnlen)sizeof(int));
      e_wsle();
      arret = true;
    }
  }
  /* -------- STARTN  SWITCH FOR STARTING VALUES OF NEWTON ITERATIONS */
  if (iwork[4] == 0) {
    startn = false;
  } else {
    startn = true;
  }
  /* -------- PARAMETER FOR DIFFERENTIAL-ALGEBRAIC COMPONENTS */
  nind1 = iwork[5];
  nind2 = iwork[6];
  nind3 = iwork[7];
  if (nind1 == 0) {
    nind1 = n;
  }
  if (nind1 + nind2 + nind3 != n) {
    s_wsle(&io___22);
    do_lio(&c__9, &c__1, " CURIOUS INPUT FOR IWORK(5,6,7)=", (ftnlen)32);
    do_lio(&c__3, &c__1, (char *)&nind1, (ftnlen)sizeof(int));
    do_lio(&c__3, &c__1, (char *)&nind2, (ftnlen)sizeof(int));
    do_lio(&c__3, &c__1, (char *)&nind3, (ftnlen)sizeof(int));
    e_wsle();
    arret = true;
  }
  /* -------- PRED   STEP SIZE CONTROL */
  if (iwork[8] <= 1) {
    pred = true;
  } else {
    pred = false;
  }
  /* -------- PARAMETER FOR SECOND ORDER EQUATIONS */
  m1 = iwork[9];
  m2 = iwork[10];
  nm1 = n - m1;
  if (m1 == 0) {
    m2 = n;
  }
  if (m2 == 0) {
    m2 = m1;
  }
  if (m1 < 0 || m2 < 0 || m1 + m2 > n) {
    s_wsle(&io___27);
    do_lio(&c__9, &c__1, " CURIOUS INPUT FOR IWORK(9,10)=", (ftnlen)31);
    do_lio(&c__3, &c__1, (char *)&m1, (ftnlen)sizeof(int));
    do_lio(&c__3, &c__1, (char *)&m2, (ftnlen)sizeof(int));
    e_wsle();
    arret = true;
  }
  /* -------- UROUND   SMALLEST NUMBER SATISFYING 1.0D0+UROUND>1.0D0 */
  if (work[1] == 0.) {
    uround = 1e-16;
  } else {
    uround = work[1];
    if (uround <= 1e-19 || uround >= 1.) {
      s_wsle(&io___29);
      do_lio(&c__9, &c__1, " COEFFICIENTS HAVE 20 DIGITS, UROUND=", (ftnlen)37);
      do_lio(&c__5, &c__1, (char *)&work[1], (ftnlen)sizeof(double));
      e_wsle();
      arret = true;
    }
  }
  /* --------- CHECK IF TOLERANCES ARE O.K. */
  if (*itol == 0) {
    if (atol[1] <= 0. || rtol[1] <= uround * 10.) {
      s_wsle(&io___30);
      do_lio(&c__9, &c__1, " TOLERANCES ARE TOO SMALL", (ftnlen)25);
      e_wsle();
      arret = true;
    }
  } else {
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      if (atol[i__] <= 0. || rtol[i__] <= uround * 10.) {
        s_wsle(&io___32);
        do_lio(&c__9, &c__1, " TOLERANCES(", (ftnlen)12);
        do_lio(&c__3, &c__1, (char *)&i__, (ftnlen)sizeof(int));
        do_lio(&c__9, &c__1, ") ARE TOO SMALL", (ftnlen)15);
        e_wsle();
        arret = true;
      }
    }
  }
  /* --------- SAFE     SAFETY FACTOR IN STEP SIZE PREDICTION */
  if (work[2] == 0.) {
    safe = .9;
  } else {
    safe = work[2];
    if (safe <= .001 || safe >= 1.) {
      s_wsle(&io___34);
      do_lio(&c__9, &c__1, " CURIOUS INPUT FOR WORK(2)=", (ftnlen)27);
      do_lio(&c__5, &c__1, (char *)&work[2], (ftnlen)sizeof(double));
      e_wsle();
      arret = true;
    }
  }
  /* ------ THET     DECIDES WHETHER THE JACOBIAN SHOULD BE RECOMPUTED; */
  if (work[3] == 0.) {
    thet = .001;
  } else {
    thet = work[3];
    if (thet >= 1.) {
      s_wsle(&io___36);
      do_lio(&c__9, &c__1, " CURIOUS INPUT FOR WORK(3)=", (ftnlen)27);
      do_lio(&c__5, &c__1, (char *)&work[3], (ftnlen)sizeof(double));
      e_wsle();
      arret = true;
    }
  }
  /* --- QUOT1 AND QUOT2: IF QUOT1 < HNEW/HOLD < QUOT2, STEP SIZE = CONST. */
  if (work[5] == 0.) {
    quot1 = 1.;
  } else {
    quot1 = work[5];
  }
  if (work[6] == 0.) {
    quot2 = 1.2;
  } else {
    quot2 = work[6];
  }
  if (quot1 > 1. || quot2 < 1.) {
    s_wsle(&io___39);
    do_lio(&c__9, &c__1, " CURIOUS INPUT FOR WORK(5,6)=", (ftnlen)29);
    do_lio(&c__5, &c__1, (char *)&quot1, (ftnlen)sizeof(double));
    do_lio(&c__5, &c__1, (char *)&quot2, (ftnlen)sizeof(double));
    e_wsle();
    arret = true;
  }
  /* -------- MAXIMAL STEP SIZE */
  if (work[7] == 0.) {
    hmax = *xend - *x;
  } else {
    hmax = work[7];
  }
  /* -------  FACL,FACR     PARAMETERS FOR STEP SIZE SELECTION */
  if (work[8] == 0.) {
    facl = 5.;
  } else {
    facl = 1. / work[8];
  }
  if (work[9] == 0.) {
    facr = .125;
  } else {
    facr = 1. / work[9];
  }
  if (facl < 1. || facr > 1.) {
    s_wsle(&io___43);
    do_lio(&c__9, &c__1, " CURIOUS INPUT WORK(8,9)=", (ftnlen)25);
    do_lio(&c__5, &c__1, (char *)&work[8], (ftnlen)sizeof(double));
    do_lio(&c__5, &c__1, (char *)&work[9], (ftnlen)sizeof(double));
    e_wsle();
    arret = true;
  }
  /* -------- PARAMETERS FOR ORDER SELECTION STRATEGY */
  if (work[10] == 0.) {
    vitu = .002;
  } else {
    vitu = work[10];
  }
  if (work[11] == 0.) {
    vitd = .8;
  } else {
    vitd = work[11];
  }
  if (work[12] == 0.) {
    hhou = 1.2;
  } else {
    hhou = work[12];
  }
  if (work[13] == 0.) {
    hhod = .8;
  } else {
    hhod = work[13];
  }
  /* *** *** *** *** *** *** *** *** *** *** *** *** *** */
  /*         COMPUTATION OF ARRAY ENTRIES */
  /* *** *** *** *** *** *** *** *** *** *** *** *** *** */
  /* ---- IMPLICIT, BANDED OR NOT ? */
  implct = *imas != 0;
  jband = *mljac < nm1;
  /* -------- COMPUTATION OF THE ROW-DIMENSIONS OF THE 2-ARRAYS --- */
  /* -- JACOBIAN  AND  MATRICES E1, E2 */
  if (jband) {
    ldjac = *mljac + *mujac + 1;
    lde1 = *mljac + ldjac;
  } else {
    *mljac = nm1;
    *mujac = nm1;
    ldjac = nm1;
    lde1 = nm1;
  }
  /* -- MASS MATRIX */
  if (implct) {
    if (*mlmas != nm1) {
      ldmas = *mlmas + *mumas + 1;
      if (jband) {
        ijob = 4;
      } else {
        ijob = 3;
      }
    } else {
      ldmas = nm1;
      ijob = 5;
    }
    /* ------ BANDWITH OF "MAS" NOT SMALLER THAN BANDWITH OF "JAC" */
    if (*mlmas > *mljac || *mumas > *mujac) {
      s_wsle(&io___54);
      do_lio(&c__9, &c__1, "BANDWITH OF \"MAS\" NOT SMALLER THAN BANDW\
ITH OF \"JAC\"",
             (ftnlen)52);
      e_wsle();
      arret = true;
    }
  } else {
    ldmas = 0;
    if (jband) {
      ijob = 2;
    } else {
      ijob = 1;
      if (n > 2 && iwork[1] != 0) {
        ijob = 7;
      }
    }
  }
  ldmas2 = max(1, ldmas);
  /* ------ HESSENBERG OPTION ONLY FOR EXPLICIT EQU. WITH FULL JACOBIAN */
  if ((implct || jband) && ijob == 7) {
    s_wsle(&io___56);
    do_lio(&c__9, &c__1, " HESSENBERG OPTION ONLY FOR EXPLICIT EQUATIONS\
 WITH FULL JACOBIAN",
           (ftnlen)65);
    e_wsle();
    arret = true;
  }
  /* ------- PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN WORK ----- */
  nns = ns * n;
  nm1ns = ns * nm1;
  nmee = (ns - 1) * nm1;
  iezz = 21;
  iey0 = iezz + nns;
  iescal = iey0 + n;
  ieff = iescal + n;
  iecon = ieff + nns;
  iejac = iecon + nns + n;
  iemas = iejac + n * ldjac;
  iee1 = iemas + nm1 * ldmas;
  iee = iee1 + nm1 * lde1;
  /* ------ TOTAL STORAGE REQUIREMENT ----------- */
  istore = iee + nmee * lde1 - 1;
  if (istore > *lwork) {
    s_wsle(&io___70);
    do_lio(&c__9, &c__1,
           " INSUFFICIENT STORAGE FOR WORK, MIN. LWORK=", (ftnlen)43);
    do_lio(&c__3, &c__1, (char *)&istore, (ftnlen)sizeof(int));
    e_wsle();
    arret = true;
  }
  /* ------- ENTRY POINTS FOR INTEGER WORKSPACE ----- */
  ieip1 = 21;
  ieip2 = ieip1 + nm1;
  ieiph = ieip2 + nm1 * (ns - 1) / 2;
  /* --------- TOTAL REQUIREMENT --------------- */
  istore = ieiph + nm1 - 1;
  if (istore > *liwork) {
    s_wsle(&io___74);
    do_lio(&c__9, &c__1,
           " INSUFF. STORAGE FOR IWORK, MIN. LIWORK=", (ftnlen)40);
    do_lio(&c__3, &c__1, (char *)&istore, (ftnlen)sizeof(int));
    e_wsle();
    arret = true;
  }
  /* ------ WHEN A FAIL HAS OCCURED, WE RETURN WITH IDID=-1 */
  if (arret) {
    *idid = -1;
    return 0;
  }
  /* -------- CALL TO CORE INTEGRATOR ------------ */
  radcov_(n, (U_fp)fcn, x, &y[1], xend, &hmax, h__, &rtol[1], &atol[1], itol,
          &nsus, (U_fp)jac, ijac, mljac, mujac, (U_fp)mas, mlmas, mumas,
          (U_fp)solout, iout, idid, &nmax, &uround, &safe, &thet, &quot1,
          &quot2, &nit, &ijob, &startn, &nind1, &nind2, &nind3, &pred, &facl,
          &facr, &m1, &m2, &nm1, &nsmin, &ns, &nns, &nm1ns, &nmee, &implct,
          &jband, &ldjac, &lde1, &ldmas2, &work[iezz], &work[iey0],
          &work[iescal], &work[ieff], &work[iejac], &work[iee1], &work[iee],
          &work[iemas], &work[iecon], &iwork[ieip1], &iwork[ieip2],
          &iwork[ieiph], &vitu, &vitd, &hhou, &hhod, &nfcn, &njac, &nstep,
          &naccpt, &nrejct, &ndec, &nsol, &rpar[1], &ipar[1]);
  iwork[13] = nsus;
  iwork[14] = nfcn;
  iwork[15] = njac;
  iwork[16] = nstep;
  iwork[17] = naccpt;
  iwork[18] = nrejct;
  iwork[19] = ndec;
  iwork[20] = nsol;
  /* ----------- RETURN ----------- */
  return 0;
} /* radau_ */

} // namespace darksun::stiff

#endif // DARKSUN_STIFF_RADAU_HPP
