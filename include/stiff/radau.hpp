//
// Created by logan on 8/10/20.
//

#ifndef STIFF_RADAU_HPP
#define STIFF_RADAU_HPP

#include "stiff/common.hpp"
#include "stiff/dc_decsol.hpp"
#include "stiff/decsol.hpp"
#include <iostream>

namespace stiff {

struct RadauWeight {
  int nn, ns;
  double xsol, hsol, c__[8];
};

struct RadauCoe3 {
  double t311, t312, t313, t321, t322, t323, t331, ti311, ti312, ti313, ti321,
      ti322, ti323, ti331, ti332, ti333;
};

struct RadauCoe5 {
  double t511, t512, t513, t514, t515, t521, t522, t523, t524, t525, t531, t532,
      t533, t534, t535, t541, t542, t543, t544, t545, t551, ti511, ti512, ti513,
      ti514, ti515, ti521, ti522, ti523, ti524, ti525, ti531, ti532, ti533,
      ti534, ti535, ti541, ti542, ti543, ti544, ti545, ti551, ti552, ti553,
      ti554, ti555;
};

struct RadauCoe7 {
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

using F_solout =
    std::function<void(int *, double *, double *, double *, double *, int *,
                       int *, int *, RadauWeight &)>;

static double c_b90 = 1.;
static double c_b100 = .8;
static double c_b140 = 9.;
static double c_b141 = .33333333333333331;

int coertv(int *nsmax, RadauCoe3 &coe3, RadauCoe5 &coe5, RadauCoe7 &coe7) {
  /* --- */
  coe3.t311 = .09123239487089294279155;
  coe3.t312 = -.141255295020954208428;
  coe3.t313 = -.03002919410514742449186;
  coe3.t321 = .2417179327071070189575;
  coe3.t322 = .204129352293799931996;
  coe3.t323 = .3829421127572619377954;
  coe3.t331 = .9660481826150929361906;
  coe3.ti311 = 4.325579890063155351024;
  coe3.ti312 = .3391992518158098695428;
  coe3.ti313 = .5417705399358748711865;
  coe3.ti321 = -4.178718591551904727346;
  coe3.ti322 = -.3276828207610623870825;
  coe3.ti323 = .4766235545005504519601;
  coe3.ti331 = -.5028726349457868759512;
  coe3.ti332 = 2.571926949855605429187;
  coe3.ti333 = -.5960392048282249249688;
  if (*nsmax <= 3) {
    return 0;
  }
  coe5.t511 = -.01251758622050104589014;
  coe5.t512 = -.01024204781790882707009;
  coe5.t513 = .04767387729029572386318;
  coe5.t514 = -.01147851525522951470794;
  coe5.t515 = -.01401985889287541028108;
  coe5.t521 = -.001491670151895382429004;
  coe5.t522 = .05017286451737105816299;
  coe5.t523 = -.09433181918161143698066;
  coe5.t524 = -.007668830749180162885157;
  coe5.t525 = .02470857842651852681253;
  coe5.t531 = -.07298187638808714862266;
  coe5.t532 = -.2305395340434179467214;
  coe5.t533 = .1027030453801258997922;
  coe5.t534 = .01939846399882895091122;
  coe5.t535 = .08180035370375117083639;
  coe5.t541 = -.3800914400035681041264;
  coe5.t542 = .3778939022488612495439;
  coe5.t543 = .4667441303324943592896;
  coe5.t544 = .4076011712801990666217;
  coe5.t545 = .1996824278868025259365;
  coe5.t551 = -.9219789736812104884883;
  coe5.ti511 = -30.04156772154440162771;
  coe5.ti512 = -13.86510785627141316518;
  coe5.ti513 = -3.480002774795185561828;
  coe5.ti514 = 1.032008797825263422771;
  coe5.ti515 = -.8043030450739899174753;
  coe5.ti521 = 5.344186437834911598895;
  coe5.ti522 = 4.593615567759161004454;
  coe5.ti523 = -3.036360323459424298646;
  coe5.ti524 = 1.05066019023145886386;
  coe5.ti525 = -.2727786118642962705386;
  coe5.ti531 = 3.748059807439804860051;
  coe5.ti532 = -3.984965736343884667252;
  coe5.ti533 = -1.044415641608018792942;
  coe5.ti534 = 1.184098568137948487231;
  coe5.ti535 = -.4499177701567803688988;
  coe5.ti541 = -33.04188021351900000806;
  coe5.ti542 = -17.37695347906356701945;
  coe5.ti543 = -.1721290632540055611515;
  coe5.ti544 = -.09916977798254264258817;
  coe5.ti545 = .5312281158383066671849;
  coe5.ti551 = -8.6114439798752919777;
  coe5.ti552 = 9.699991409528808231336;
  coe5.ti553 = 1.914728639696874284851;
  coe5.ti554 = 2.418692006084940026427;
  coe5.ti555 = -1.047463487935337418694;
  if (*nsmax <= 5) {
    return 0;
  }
  coe7.t711 = -.002153754627310526422828;
  coe7.t712 = .02156755135132077338691;
  coe7.t713 = .008783567925144144407326;
  coe7.t714 = -.004055161452331023898198;
  coe7.t715 = .004427232753268285479678;
  coe7.t716 = -.001238646187952874056377;
  coe7.t717 = -.002760617480543852499548;
  coe7.t721 = .001600025077880428526831;
  coe7.t722 = -.03813164813441154669442;
  coe7.t723 = -.02152556059400687552385;
  coe7.t724 = .008415568276559589237177;
  coe7.t725 = -.004031949570224549492304;
  coe7.t726 = -6.666635339396338181761e-5;
  coe7.t727 = .003185474825166209848748;
  coe7.t731 = -.00405910730194768309165;
  coe7.t732 = .05739650893938171539757;
  coe7.t733 = .05885052920842679105612;
  coe7.t734 = -.008560431061603432060177;
  coe7.t735 = -.006923212665023908924141;
  coe7.t736 = -.002352180982943338340535;
  coe7.t737 = 4.169077725297562691409e-4;
  coe7.t741 = -.01575048807937684420346;
  coe7.t742 = -.03821469359696835048464;
  coe7.t743 = -.1657368112729438512412;
  coe7.t744 = -.03737124230238445741907;
  coe7.t745 = .008239007298507719404499;
  coe7.t746 = .003115071152346175252726;
  coe7.t747 = .02511660491343882192836;
  coe7.t751 = -.1129776610242208076086;
  coe7.t752 = -.2491742124652636863308;
  coe7.t753 = .2735633057986623212132;
  coe7.t754 = .005366761379181770094279;
  coe7.t755 = .1932111161012620144312;
  coe7.t756 = .1017177324817151468081;
  coe7.t757 = .09504502035604622821039;
  coe7.t761 = -.4583810431839315010281;
  coe7.t762 = .5315846490836284292051;
  coe7.t763 = .4863228366175728940567;
  coe7.t764 = .5265742264584492629141;
  coe7.t765 = .2755343949896258141929;
  coe7.t766 = .5217519452747652852946;
  coe7.t767 = .1280719446355438944141;
  coe7.t771 = -.8813915783538183763135;
  coe7.ti711 = -258.1319263199822292761;
  coe7.ti712 = -189.073763081398508952;
  coe7.ti713 = -49.08731481793013119445;
  coe7.ti714 = -4.110647469661428418112;
  coe7.ti715 = -4.053447889315563304175;
  coe7.ti716 = 3.112755366607346076554;
  coe7.ti717 = -1.646774913558444650169;
  coe7.ti721 = -3.007390169451292131731;
  coe7.ti722 = -11.01586607876577132911;
  coe7.ti723 = 1.487799456131656281486;
  coe7.ti724 = 2.130388159559282459432;
  coe7.ti725 = -1.816141086817565624822;
  coe7.ti726 = 1.134325587895161100083;
  coe7.ti727 = -.414699045943303531993;
  coe7.ti731 = -8.441963188321084681757;
  coe7.ti732 = -.6505252740575150028169;
  coe7.ti733 = 6.940670730369876478804;
  coe7.ti734 = -3.205047525597898431565;
  coe7.ti735 = 1.071280943546478589783;
  coe7.ti736 = -.354850749121622187973;
  coe7.ti737 = .09198549132786554154409;
  coe7.ti741 = 74.67833223502269977153;
  coe7.ti742 = 87.40858897990081640204;
  coe7.ti743 = 4.024158737379997877014;
  coe7.ti744 = -3.714806315158364186639;
  coe7.ti745 = -3.430093985982317350741;
  coe7.ti746 = 2.696604809765312378853;
  coe7.ti747 = -.9386927436075461933568;
  coe7.ti751 = 58.35652885190657724237;
  coe7.ti752 = -10.06877395780018096325;
  coe7.ti753 = -30.36638884256667120811;
  coe7.ti754 = -1.020020865184865985027;
  coe7.ti755 = -.1124175003784249621267;
  coe7.ti756 = 1.8906408310003776228;
  coe7.ti757 = -.9716486393831482282172;
  coe7.ti761 = -299.1862480282520966786;
  coe7.ti762 = -243.0407453687447911819;
  coe7.ti763 = -48.77710407803786921219;
  coe7.ti764 = -2.03867190574193440528;
  coe7.ti765 = 1.673560239861084944268;
  coe7.ti766 = -1.087374032057106164456;
  coe7.ti767 = .9019382492960993738427;
  coe7.ti771 = -93.07650289743530591157;
  coe7.ti772 = 23.88163105628114427703;
  coe7.ti773 = 39.2788807308138438271;
  coe7.ti774 = 14.38891568549108006988;
  coe7.ti775 = -3.510438399399361221087;
  coe7.ti776 = 4.863284885566180701215;
  coe7.ti777 = -2.2464827295912399164;
  return 0;
}

int coercv(int *ns, double *c__, double *dd, double *u1, double *alph,
           double *beta) {
  /* System generated locals */
  double d__1, d__2;

  /* Local variables */
  double sq6, st9, bet, alp, cno;

  /* Parameter adjustments */
  --beta;
  --alph;
  --dd;

  /* Function Body */
  c__[0] = 0.;
  c__[*ns] = 1.;
  switch (*ns) {
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
  }
L11:
  return 0;
L1:
  c__[1] = 1.;
  *u1 = 1.;
  dd[1] = -1.;
  return 0;
L3:
  sq6 = sqrt(6.);
  c__[1] = (4. - sq6) / 10.;
  c__[2] = (sq6 + 4.) / 10.;
  st9 = pow(c_b140, c_b141);
  *u1 = (st9 * (st9 - 1) + 6.) / 30.;
  alp = (12. - st9 * (st9 - 1)) / 60.;
  bet = st9 * (st9 + 1) * sqrt(3.) / 60.;
  /* Computing 2nd power */
  d__1 = alp;
  /* Computing 2nd power */
  d__2 = bet;
  cno = d__1 * d__1 + d__2 * d__2;
  *u1 = 1. / *u1;
  alph[1] = alp / cno;
  beta[1] = bet / cno;
  return 0;
L5:
  c__[1] = .05710419611451768219312;
  c__[2] = .27684301363812382768;
  c__[3] = .5835904323689168200567;
  c__[4] = .8602401356562194478479;
  dd[1] = -27.78093394406463730479;
  dd[2] = 3.641478498049213152712;
  dd[3] = -1.252547721169118720491;
  dd[4] = .5920031671845428725662;
  dd[5] = -.2;
  *u1 = 6.286704751729276645173;
  alph[1] = 3.655694325463572258243;
  beta[1] = 6.543736899360077294021;
  alph[2] = 5.70095329867178941917;
  beta[2] = 3.210265600308549888425;
  return 0;
L7:
  c__[1] = .02931642715978489197205;
  c__[2] = .14807859966848429185;
  c__[3] = .3369846902811542990971;
  c__[4] = .5586715187715501320814;
  c__[5] = .7692338620300545009169;
  c__[6] = .9269456713197411148519;
  dd[1] = -54.37443689412861451458;
  dd[2] = 7.000024004259186512041;
  dd[3] = -2.355661091987557192256;
  dd[4] = 1.132289066106134386384;
  dd[5] = -.6468913267673587118673;
  dd[6] = .3875333853753523774248;
  dd[7] = -.1428571428571428571429;
  *u1 = 8.936832788405216337302;
  alph[1] = 4.378693561506806002523;
  beta[1] = 10.16969328379501162732;
  alph[2] = 7.141055219187640105775;
  beta[2] = 6.623045922639275970621;
  alph[3] = 8.511834825102945723051;
  beta[3] = 3.281013624325058830036;
  return 0;
}

double contra(int *i__, double *x, double *cont, int * /*lrc*/,
              const RadauWeight &weight) {
  /* System generated locals */
  double ret_val;

  /* Local variables */
  int k;
  double s;

  /* ---------------------------------------------------------- */
  /*     THIS FUNCTION CAN BE USED FOR CONINUOUS OUTPUT. IT PROVIDES AN */
  /*     APPROXIMATION TO THE I-TH COMPONENT OF THE SOLUTION AT X. */
  /*     IT GIVES THE VALUE OF THE COLLOCATION POLYNOMIAL, DEFINED FOR */
  /*     THE LAST SUCCESSFULLY COMPUTED STEP (BY RADAU). */
  /* ---------------------------------------------------------- */
  /* Parameter adjustments */
  --cont;

  /* Function Body */
  s = (*x - weight.xsol) / weight.hsol + 1.;
  ret_val = cont[*i__ + weight.ns * weight.nn];
  for (k = weight.ns - 1; k >= 0; --k) {
    ret_val =
        cont[*i__ + k * weight.nn] + (s - weight.c__[weight.ns - k]) * ret_val;
  }
  return ret_val;
}

//===========================================================================
//---- Radau Core Integrator ------------------------------------------------
//===========================================================================

int radcov(int *n, F_fcn fcn, double *x, double *y, double *xend, double *hmax,
           double *h__, double *rtol, double *atol, int *itol, int *ns,
           F_jac jac, int *ijac, int *mljac, int *mujac, F_mas mas, int *mlmas,
           int *mumas, F_solout solout, int *iout, int *idid, int *nmax,
           double *uround, double *safe, double *thet, double *quot1,
           double *quot2, int *nit1, int *ijob, bool *startn, int *nind1,
           int *nind2, int *nind3, bool *pred, double *facl, double *facr,
           int *m1, int *m2, int *nm1, int *nsmin, int *nsmax, int * /*nnms*/,
           int * /*nm1ns*/, int * /*nmee*/, bool *implct, bool *banded,
           int *ldjac, int *lde1, int *ldmas, double *zz, double *y0,
           double *scal, double *ff, double *fjac, double *e1, double *ee2,
           double *fmas, double *cont, int *ip1, int *ip2, int *iphes,
           double *vitu, double *vitd, double *hhou, double *hhod, int *nfcn,
           int *njac, int *nstep, int *naccpt, int *nrejct, int *ndec,
           int *nsol) {

  /* System generated locals */
  int fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e1_dim1, e1_offset,
      ee2_dim1, ee2_offset, i__1, i__2, i__3, i__4;
  double d__1, d__2, d__3, d__4, d__5;

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
  /* Parameter adjustments */
  --scal;
  --y0;
  --y;
  --rtol;
  --atol;
  --iphes;
  --ip1;
  --cont;
  --ff;
  --zz;
  --ip2;
  fjac_dim1 = *ldjac;
  fjac_offset = 1 + fjac_dim1;
  fjac -= fjac_offset;
  ee2_dim1 = *lde1;
  ee2_offset = 1 + ee2_dim1;
  ee2 -= ee2_offset;
  e1_dim1 = *lde1;
  e1_offset = 1 + e1_dim1;
  e1 -= e1_offset;
  fmas_dim1 = *ldmas;
  fmas_offset = 1 + fmas_dim1;
  fmas -= fmas_offset;

  /* Function Body */
  index1 = *nind1 != 0;
  index2 = *nind2 != 0;
  index3 = *nind3 != 0;
  /* ------- COMPUTE MASS MATRIX FOR IMPLICIT CASE ---------- */
  if (*implct) {
    mas(nm1, &fmas[fmas_offset], ldmas);
  }
  variab = *nsmin < *nsmax;
  /* ---------- CONSTANTS --------- */
  expo = 1. / (*ns + 1.);
  sq6 = sqrt(6.0);
  c31 = (4.0 - sq6) / 10.;
  c32 = (sq6 + 4.) / 10.;
  c31m1 = c31 - 1.;
  c32m1 = c32 - 1.;
  c31mc2 = c31 - c32;
  dd1 = -(sq6 * 7. + 13.) / 3.;
  dd2 = (sq6 * 7. - 13.) / 3.;
  dd3 = -.33333333333333331;
  n2 = *n << 1;
  n3 = *n * 3;
  n4 = *n << 2;
  n5 = *n * 5;
  n6 = *n * 6;
  unexp = false;
  unexn = false;
  change = false;
  ikeep = 0;
  ichan = 0;
  theta = 0.;
  thetat = 0.;

  RadauWeight weight_1;
  weight_1.nn = *n;
  nns = *n * *ns;
  weight_1.ns = *ns;
  lrc = nns + *n;

  RadauCoe3 coe3;
  RadauCoe5 coe5;
  RadauCoe7 coe7;
  coertv(nsmax, coe3, coe5, coe7);
  coercv(ns, weight_1.c__, dd, &u1, alph, beta);
  if (*m1 > 0) {
    *ijob += 10;
  }
  d__1 = *xend - *x;
  posneg = std::copysign(c_b90, d__1);
  /* Computing std::min */
  d__2 = std::abs(*hmax), d__3 = (d__1 = *xend - *x, std::abs(d__1));
  hmaxn = std::min(d__2, d__3);
  if (std::abs(*h__) <= *uround * 10.) {
    *h__ = 1e-6;
  }
  /* Computing std::min */
  d__1 = std::abs(*h__);
  *h__ = std::min(d__1, hmaxn);
  *h__ = std::copysign(*h__, posneg);
  hold = *h__;
  reject = false;
  first = true;
  last = false;
  if ((*x + *h__ * 1.0001 - *xend) * posneg >= 0.) {
    *h__ = *xend - *x;
    last = true;
  }
  hopt = *h__;
  faccon = 1.;
  nsing = 0;
  xold = *x;
  if (*iout != 0) {
    irtrn = 1;
    nrsol = 1;
    xosol = xold;
    weight_1.xsol = *x;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      cont[i__] = y[i__];
    }
    nsolu = *n;
    weight_1.hsol = hold;
    solout(&nrsol, &xosol, &weight_1.xsol, &y[1], &cont[1], &lrc, &nsolu,
           &irtrn, weight_1);
    if (irtrn < 0) {
      goto L179;
    }
  }
  StiffLinAlg linal_1;
  linal_1.mle = *mljac;
  linal_1.mue = *mujac;
  linal_1.mbjac = *mljac + *mujac + 1;
  linal_1.mbb = *mlmas + *mumas + 1;
  linal_1.mdiag = linal_1.mle + linal_1.mue + 1;
  linal_1.mdiff = linal_1.mle + linal_1.mue - *mumas;
  linal_1.mbdiag = *mumas + 1;
  expmns = (*ns + 1.) / (*ns * 2.);
  quott = atol[1] / rtol[1];
  if (*itol == 0) {
    rtol1 = pow(rtol[1], expmns) * .1;
    atol1 = rtol1 * quott;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      scal[i__] = atol1 + rtol1 * (d__1 = y[i__], std::abs(d__1));
    }
  } else {
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      quott = atol[i__] / rtol[i__];
      rtol1 = pow(rtol[i__], expmns) * .1;
      atol1 = rtol1 * quott;
      scal[i__] = atol1 + rtol1 * (d__1 = y[i__], std::abs(d__1));
    }
  }
  hhfac = *h__;
  fcn(n, x, &y[1], &y0[1]);
  ++(*nfcn);
  /* --- BASIC INTEGRATION STEP */
L10:
  /* *** *** *** *** *** *** *** */
  /*  COMPUTATION OF THE JACOBIAN */
  /* *** *** *** *** *** *** *** */
  ++(*njac);
  if (*ijac == 0) {
    /* --- COMPUTE JACOBIAN MATRIX NUMERICALLY */
    if (*banded) {
      /* --- JACOBIAN IS BANDED */
      mujacp = *mujac + 1;
      md = std::min(linal_1.mbjac, *m2);
      i__1 = *m1 / *m2 + 1;
      for (mm = 1; mm <= i__1; ++mm) {
        i__2 = md;
        for (k = 1; k <= i__2; ++k) {
          j = k + (mm - 1) * *m2;
        L12:
          ff[j] = y[j];
          /* Computing std::max */
          d__2 = 1e-5, d__3 = (d__1 = y[j], std::abs(d__1));
          ff[j + *n] = sqrt(*uround * std::max(d__2, d__3));
          y[j] += ff[j + *n];
          j += md;
          if (j <= mm * *m2) {
            goto L12;
          }
          fcn(n, x, &y[1], &cont[1]);
          j = k + (mm - 1) * *m2;
          j1 = k;
          /* Computing std::max */
          i__3 = 1, i__4 = j1 - *mujac;
          lbeg = std::max(i__3, i__4) + *m1;
        L14:
          /* Computing std::min */
          i__3 = *m2, i__4 = j1 + *mljac;
          lend = std::min(i__3, i__4) + *m1;
          y[j] = ff[j];
          mujacj = mujacp - j1 - *m1;
          i__3 = lend;
          for (l = lbeg; l <= i__3; ++l) {
            fjac[l + mujacj + j * fjac_dim1] = (cont[l] - y0[l]) / ff[j + *n];
          }
          j += md;
          j1 += md;
          lbeg = lend + 1;
          if (j <= mm * *m2) {
            goto L14;
          }
        }
      }
    } else {
      /* --- JACOBIAN IS FULL */
      i__1 = *n;
      for (i__ = 1; i__ <= i__1; ++i__) {
        ysafe = y[i__];
        /* Computing std::max */
        d__1 = 1e-5, d__2 = std::abs(ysafe);
        delt = sqrt(*uround * std::max(d__1, d__2));
        y[i__] = ysafe + delt;
        fcn(n, x, &y[1], &cont[1]);
        i__2 = *n;
        for (j = *m1 + 1; j <= i__2; ++j) {
          fjac[j - *m1 + i__ * fjac_dim1] = (cont[j] - y0[j]) / delt;
        }
        y[i__] = ysafe;
      }
    }
  } else {
    /* --- COMPUTE JACOBIAN MATRIX ANALYTICALLY */
    jac(n, x, &y[1], &fjac[fjac_offset], ldjac);
  }
  caljac = true;
  calhes = true;
L20:
  /* --- CHANGE THE ORDER HERE IF NECESSARY */
  if (variab) {
    nsnew = *ns;
    ++ichan;
    hquot = *h__ / hold;
    /* Computing std::min */
    /* Computing std::max */
    d__3 = theta, d__4 = thetat * 0.5;
    d__1 = 10., d__2 = std::max(d__3, d__4);
    thetat = std::min(d__1, d__2);
    if (newt > 1 && thetat <= *vitu && hquot < *hhou && hquot > *hhod) {
      /* Computing std::min */
      i__1 = *nsmax, i__2 = *ns + 2;
      nsnew = std::min(i__1, i__2);
    }
    if (thetat >= *vitd || unexp) {
      /* Computing std::max */
      i__1 = *nsmin, i__2 = *ns - 2;
      nsnew = std::max(i__1, i__2);
    }
    if (ichan >= 1 && unexn) {
      /* Computing std::max */
      i__1 = *nsmin, i__2 = *ns - 2;
      nsnew = std::max(i__1, i__2);
    }
    if (ichan <= 10) {
      nsnew = std::min(*ns, nsnew);
    }
    change = *ns != nsnew;
    unexn = false;
    unexp = false;
    if (change) {
      *ns = nsnew;
      ichan = 1;
      nns = *n * *ns;
      weight_1.ns = *ns;
      lrc = nns + *n;
      coercv(ns, weight_1.c__, dd, &u1, alph, beta);
      expo = 1. / (*ns + 1.);
      expmns = (*ns + 1.) / (*ns * 2.);
      rtol1 = pow(rtol[1], expmns) * .1;
      atol1 = rtol1 * quott;
      if (*itol == 0) {
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
          scal[i__] = atol1 + rtol1 * (d__1 = y[i__], std::abs(d__1));
        }
      } else {
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
          quott = atol[i__] / rtol[i__];
          rtol1 = pow(rtol[i__], expmns) * .1;
          atol1 = rtol1 * quott;
          scal[i__] = atol1 + rtol1 * (d__1 = y[i__], std::abs(d__1));
        }
      }
    }
  }
  /* --- COMPUTE THE MATRICES E1 AND E2 AND THEIR DECOMPOSITIONS */
  fac1 = u1 / *h__;
  decomr(n, &fjac[fjac_offset], ldjac, &fmas[fmas_offset], ldmas, mlmas, mumas,
         m1, m2, nm1, &fac1, &e1[e1_offset], lde1, &ip1[1], &ier, ijob, &calhes,
         &iphes[1], linal_1);
  if (ier != 0) {
    goto L78;
  }
  i__1 = (*ns - 1) / 2;
  for (k = 1; k <= i__1; ++k) {
    alphn[k - 1] = alph[k - 1] / *h__;
    betan[k - 1] = beta[k - 1] / *h__;
    // v-- this was: iad = (k - 1 << 1) * *nm1 + 1;
    iad = ((k - 1) << 1) * *nm1 + 1;
    decomc(n, &fjac[fjac_offset], ldjac, &fmas[fmas_offset], ldmas, mlmas,
           mumas, m1, m2, nm1, &alphn[k - 1], &betan[k - 1],
           &ee2[iad * ee2_dim1 + 1], &ee2[(iad + *nm1) * ee2_dim1 + 1], lde1,
           &ip2[(k - 1) * *nm1 + 1], &ier, ijob, linal_1);
    if (ier != 0) {
      goto L78;
    }
  }
  ++(*ndec);
L30:
  if (variab && ikeep == 1) {
    ++ichan;
    ikeep = 0;
    if (ichan >= 10 && *ns < *nsmax) {
      goto L20;
    }
  }
  ++(*nstep);
  if (*nstep > *nmax) {
    goto L178;
  }
  if (std::abs(*h__) * .1 <= std::abs(*x) * *uround) {
    goto L177;
  }
  if (index2) {
    i__1 = *nind1 + *nind2;
    for (i__ = *nind1 + 1; i__ <= i__1; ++i__) {
      scal[i__] /= hhfac;
    }
  }
  if (index3) {
    i__1 = *nind1 + *nind2 + *nind3;
    for (i__ = *nind1 + *nind2 + 1; i__ <= i__1; ++i__) {
      scal[i__] /= hhfac * hhfac;
    }
  }
  xph = *x + *h__;
  /* *** *** *** *** *** *** *** */
  /* *** *** *** *** *** *** *** */
  if (*ns == 3) {
    /* *** *** *** *** *** *** *** */
    /* *** *** *** *** *** *** *** */
    if (first || *startn || change) {
      i__1 = nns;
      for (i__ = 1; i__ <= i__1; ++i__) {
        zz[i__] = 0.;
        ff[i__] = 0.;
      }
    } else {
      hquot = *h__ / hold;
      c3q = hquot;
      c1q = c31 * c3q;
      c2q = c32 * c3q;
      i__1 = *n;
      for (i__ = 1; i__ <= i__1; ++i__) {
        ak1 = cont[i__ + *n];
        ak2 = cont[i__ + n2];
        ak3 = cont[i__ + n3];
        z1i = c1q * (ak1 + (c1q - c32m1) * (ak2 + (c1q - c31m1) * ak3));
        z2i = c2q * (ak1 + (c2q - c32m1) * (ak2 + (c2q - c31m1) * ak3));
        z3i = c3q * (ak1 + (c3q - c32m1) * (ak2 + (c3q - c31m1) * ak3));
        zz[i__] = z1i;
        zz[i__ + *n] = z2i;
        zz[i__ + n2] = z3i;
        ff[i__] = coe3.ti311 * z1i + coe3.ti312 * z2i + coe3.ti313 * z3i;
        ff[i__ + *n] = coe3.ti321 * z1i + coe3.ti322 * z2i + coe3.ti323 * z3i;
        ff[i__ + n2] = coe3.ti331 * z1i + coe3.ti332 * z2i + coe3.ti333 * z3i;
      }
    }
    /* *** *** *** *** *** *** *** */
    /*  LOOP FOR THE SIMPLIFIED NEWTON ITERATION */
    /* *** *** *** *** *** *** *** */
    newt = 0;
    nit = *nit1;
    expmi = 1. / expmns;
    /* Computing std::max */
    /* Computing std::min */
    d__5 = expmi - 1.;
    d__3 = 0.03, d__4 = pow(rtol1, d__5);
    d__1 = *uround * 10 / rtol1, d__2 = std::min(d__3, d__4);
    fnewt = std::max(d__1, d__2);
    d__1 = std::max(faccon, *uround);
    faccon = pow(d__1, c_b100);
    theta = std::abs(*thet);
  L40:
    if (newt >= nit) {
      goto L78;
    }
    /* ---     COMPUTE THE RIGHT-HAND SIDE */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      cont[i__] = y[i__] + zz[i__];
    }
    d__1 = *x + c31 * *h__;
    fcn(n, &d__1, &cont[1], &zz[1]);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      cont[i__] = y[i__] + zz[i__ + *n];
    }
    d__1 = *x + c32 * *h__;
    fcn(n, &d__1, &cont[1], &zz[*n + 1]);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      cont[i__] = y[i__] + zz[i__ + n2];
    }
    fcn(n, &xph, &cont[1], &zz[n2 + 1]);
    *nfcn += 3;
    /* ---     SOLVE THE LINEAR SYSTEMS */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      a1 = zz[i__];
      a2 = zz[i__ + *n];
      a3 = zz[i__ + n2];
      zz[i__] = coe3.ti311 * a1 + coe3.ti312 * a2 + coe3.ti313 * a3;
      zz[i__ + *n] = coe3.ti321 * a1 + coe3.ti322 * a2 + coe3.ti323 * a3;
      zz[i__ + n2] = coe3.ti331 * a1 + coe3.ti332 * a2 + coe3.ti333 * a3;
    }
    slvrad(n, &fjac[fjac_offset], ldjac, mljac, mujac, &fmas[fmas_offset],
           ldmas, mlmas, mumas, m1, m2, nm1, &fac1, alphn, betan,
           &e1[e1_offset], &ee2[ee2_offset], &ee2[(*nm1 + 1) * ee2_dim1 + 1],
           lde1, &zz[1], &zz[*n + 1], &zz[n2 + 1], &ff[1], &ff[*n + 1],
           &ff[n2 + 1], &cont[1], &ip1[1], &ip2[1], &iphes[1], &ier, ijob,
           linal_1);
    ++(*nsol);
    ++newt;
    dyno = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      denom = scal[i__];
      /* Computing 2nd power */
      d__1 = zz[i__] / denom;
      /* Computing 2nd power */
      d__2 = zz[i__ + *n] / denom;
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
        dyth = faccon * dyno * pow(theta, i__1) / fnewt;
        if (dyth >= 1.) {
          /* Computing std::max */
          d__1 = 1e-4, d__2 = std::min(20., dyth);
          qnewt = std::max(d__1, d__2);
          d__1 = -1. / (nit + 4. - 1 - newt);
          hhfac = pow(qnewt, d__1) * .8;
          *h__ = hhfac * *h__;
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
    dynold = std::max(dyno, *uround);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      in = i__ + *n;
      in2 = i__ + n2;
      f1i = ff[i__] + zz[i__];
      f2i = ff[in] + zz[in];
      f3i = ff[in2] + zz[in2];
      ff[i__] = f1i;
      ff[in] = f2i;
      ff[in2] = f3i;
      zz[i__] = coe3.t311 * f1i + coe3.t312 * f2i + coe3.t313 * f3i;
      zz[in] = coe3.t321 * f1i + coe3.t322 * f2i + coe3.t323 * f3i;
      zz[in2] = coe3.t331 * f1i + f2i;
    }
    if (faccon * dyno > fnewt) {
      goto L40;
    }
    /* --- ERROR ESTIMATION */
    estrad(n, &fjac[fjac_offset], ldjac, mljac, mujac, &fmas[fmas_offset],
           ldmas, mlmas, mumas, h__, &dd1, &dd2, &dd3, fcn, nfcn, &y0[1], &y[1],
           ijob, x, m1, m2, nm1, &e1[e1_offset], lde1, &zz[1], &zz[*n + 1],
           &zz[n2 + 1], &cont[1], &ff[1], &ff[*n + 1], &ip1[1], &iphes[1],
           &scal[1], &err, &first, &reject, &fac1, linal_1);
    /*       --- COMPUTE FINITE DIFFERENCES FOR DENSE OUTPUT */
    if (err < 1.) {
      i__1 = *n;
      for (i__ = 1; i__ <= i__1; ++i__) {
        y[i__] += zz[i__ + n2];
        z2i = zz[i__ + *n];
        z1i = zz[i__];
        cont[i__ + *n] = (z2i - zz[i__ + n2]) / c32m1;
        ak = (z1i - z2i) / c31mc2;
        acont3 = z1i / c31;
        acont3 = (ak - acont3) / c32;
        cont[i__ + n2] = (ak - cont[i__ + *n]) / c31m1;
        cont[i__ + n3] = cont[i__ + n2] - acont3;
      }
    }
    /* *** *** *** *** *** *** *** */
    /* *** *** *** *** *** *** *** */
  } else {
    if (*ns == 5) {
      /* *** *** *** *** *** *** *** */
      /* *** *** *** *** *** *** *** */
      if (first || *startn || change) {
        i__1 = nns;
        for (i__ = 1; i__ <= i__1; ++i__) {
          zz[i__] = 0.;
          ff[i__] = 0.;
        }
      } else {
        hquot = *h__ / hold;
        i__1 = *ns;
        for (k = 1; k <= i__1; ++k) {
          ccq = weight_1.c__[k] * hquot;
          i__2 = *n;
          for (i__ = 1; i__ <= i__2; ++i__) {
            val = cont[i__ + *ns * *n];
            for (l = *ns - 1; l >= 1; --l) {
              val =
                  cont[i__ + l * *n] + (ccq - weight_1.c__[*ns - l] + 1.) * val;
            }
            zz[i__ + (k - 1) * *n] = ccq * val;
          }
        }
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
          z1i = zz[i__];
          z2i = zz[i__ + *n];
          z3i = zz[i__ + n2];
          z4i = zz[i__ + n3];
          z5i = zz[i__ + n4];
          ff[i__] = coe5.ti511 * z1i + coe5.ti512 * z2i + coe5.ti513 * z3i +
                    coe5.ti514 * z4i + coe5.ti515 * z5i;
          ff[i__ + *n] = coe5.ti521 * z1i + coe5.ti522 * z2i +
                         coe5.ti523 * z3i + coe5.ti524 * z4i + coe5.ti525 * z5i;
          ff[i__ + n2] = coe5.ti531 * z1i + coe5.ti532 * z2i +
                         coe5.ti533 * z3i + coe5.ti534 * z4i + coe5.ti535 * z5i;
          ff[i__ + n3] = coe5.ti541 * z1i + coe5.ti542 * z2i +
                         coe5.ti543 * z3i + coe5.ti544 * z4i + coe5.ti545 * z5i;
          ff[i__ + n4] = coe5.ti551 * z1i + coe5.ti552 * z2i +
                         coe5.ti553 * z3i + coe5.ti554 * z4i + coe5.ti555 * z5i;
        }
      }
      /* *** *** *** *** *** *** *** */
      /*  LOOP FOR THE SIMPLIFIED NEWTON ITERATION */
      /* *** *** *** *** *** *** *** */
      newt = 0;
      nit = *nit1 + 5;
      expmi = 1. / expmns;
      /* Computing std::max */
      /* Computing std::min */
      d__5 = expmi - 1.;
      d__3 = .03, d__4 = pow(rtol1, d__5);
      d__1 = *uround * 10 / rtol1, d__2 = std::min(d__3, d__4);
      fnewt = std::max(d__1, d__2);
      d__1 = std::max(faccon, *uround);
      faccon = pow(d__1, c_b100);
      theta = std::abs(*thet);
    L140:
      if (newt >= nit) {
        goto L78;
      }
      /* ---     COMPUTE THE RIGHT-HAND SIDE */
      i__1 = *ns - 1;
      for (k = 0; k <= i__1; ++k) {
        iadd = k * *n;
        i__2 = *n;
        for (i__ = 1; i__ <= i__2; ++i__) {
          cont[i__] = y[i__] + zz[iadd + i__];
        }
        d__1 = *x + weight_1.c__[k + 1] * *h__;
        fcn(n, &d__1, &cont[1], &zz[iadd + 1]);
      }
      *nfcn += *ns;
      /* ---     SOLVE THE LINEAR SYSTEMS */
      i__1 = *n;
      for (i__ = 1; i__ <= i__1; ++i__) {
        z1i = zz[i__];
        z2i = zz[i__ + *n];
        z3i = zz[i__ + n2];
        z4i = zz[i__ + n3];
        z5i = zz[i__ + n4];
        zz[i__] = coe5.ti511 * z1i + coe5.ti512 * z2i + coe5.ti513 * z3i +
                  coe5.ti514 * z4i + coe5.ti515 * z5i;
        zz[i__ + *n] = coe5.ti521 * z1i + coe5.ti522 * z2i + coe5.ti523 * z3i +
                       coe5.ti524 * z4i + coe5.ti525 * z5i;
        zz[i__ + n2] = coe5.ti531 * z1i + coe5.ti532 * z2i + coe5.ti533 * z3i +
                       coe5.ti534 * z4i + coe5.ti535 * z5i;
        zz[i__ + n3] = coe5.ti541 * z1i + coe5.ti542 * z2i + coe5.ti543 * z3i +
                       coe5.ti544 * z4i + coe5.ti545 * z5i;
        zz[i__ + n4] = coe5.ti551 * z1i + coe5.ti552 * z2i + coe5.ti553 * z3i +
                       coe5.ti554 * z4i + coe5.ti555 * z5i;
      }
      slvrar(n, &fjac[fjac_offset], ldjac, mljac, mujac, &fmas[fmas_offset],
             ldmas, mlmas, mumas, m1, m2, nm1, &fac1, &e1[e1_offset], lde1,
             &zz[1], &ff[1], &ip1[1], &iphes[1], &ier, ijob, linal_1);
      for (k = 1; k <= 2; ++k) {
        // v---- this was: iad = (k - 1 << 1) * *nm1 + 1;
        iad = ((k - 1) << 1) * *nm1 + 1;
        slvrai(n, &fjac[fjac_offset], ldjac, mljac, mujac, &fmas[fmas_offset],
               ldmas, mlmas, mumas, m1, m2, nm1, &alphn[k - 1], &betan[k - 1],
               &ee2[iad * ee2_dim1 + 1], &ee2[(iad + *nm1) * ee2_dim1 + 1],
               lde1, &zz[((k << 1) - 1) * *n + 1], &zz[(k << 1) * *n + 1],
               &ff[((k << 1) - 1) * *n + 1], &ff[(k << 1) * *n + 1], &cont[1],
               &ip2[(k - 1) * *nm1 + 1], &iphes[1], &ier, ijob, linal_1);
      }
      ++(*nsol);
      ++newt;
      dyno = 0.;
      i__1 = *n;
      for (i__ = 1; i__ <= i__1; ++i__) {
        denom = scal[i__];
        i__2 = *ns - 1;
        for (k = 0; k <= i__2; ++k) {
          /* Computing 2nd power */
          d__1 = zz[i__ + k * *n] / denom;
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
          dyth = faccon * dyno * pow(theta, i__1) / fnewt;
          if (dyth >= 1.) {
            /* Computing std::max */
            d__1 = 1e-4, d__2 = std::min(20., dyth);
            qnewt = std::max(d__1, d__2);
            d__1 = -1. / (nit + 4. - 1 - newt);
            hhfac = pow(qnewt, d__1) * .8;
            *h__ = hhfac * *h__;
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
      dynold = std::max(dyno, *uround);
      i__1 = *n;
      for (i__ = 1; i__ <= i__1; ++i__) {
        z1i = ff[i__] + zz[i__];
        z2i = ff[i__ + *n] + zz[i__ + *n];
        z3i = ff[i__ + n2] + zz[i__ + n2];
        z4i = ff[i__ + n3] + zz[i__ + n3];
        z5i = ff[i__ + n4] + zz[i__ + n4];
        ff[i__] = z1i;
        ff[i__ + *n] = z2i;
        ff[i__ + n2] = z3i;
        ff[i__ + n3] = z4i;
        ff[i__ + n4] = z5i;
        zz[i__] = coe5.t511 * z1i + coe5.t512 * z2i + coe5.t513 * z3i +
                  coe5.t514 * z4i + coe5.t515 * z5i;
        zz[i__ + *n] = coe5.t521 * z1i + coe5.t522 * z2i + coe5.t523 * z3i +
                       coe5.t524 * z4i + coe5.t525 * z5i;
        zz[i__ + n2] = coe5.t531 * z1i + coe5.t532 * z2i + coe5.t533 * z3i +
                       coe5.t534 * z4i + coe5.t535 * z5i;
        zz[i__ + n3] = coe5.t541 * z1i + coe5.t542 * z2i + coe5.t543 * z3i +
                       coe5.t544 * z4i + coe5.t545 * z5i;
        zz[i__ + n4] = coe5.t551 * z1i + z2i + z4i;
      }
      if (faccon * dyno > fnewt) {
        goto L140;
      }
      /* --- ERROR ESTIMATION */
      estrav(n, &fjac[fjac_offset], ldjac, mljac, mujac, &fmas[fmas_offset],
             ldmas, mlmas, mumas, h__, dd, fcn, nfcn, &y0[1], &y[1], ijob, x,
             m1, m2, nm1, ns, &nns, &e1[e1_offset], lde1, &zz[1], &cont[1],
             &ff[1], &ip1[1], &iphes[1], &scal[1], &err, &first, &reject, &fac1,
             linal_1);
      /*       --- COMPUTE FINITE DIFFERENCES FOR DENSE OUTPUT */
      if (err < 1.) {
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
          y[i__] += zz[i__ + n4];
          cont[i__ + n5] = zz[i__] / weight_1.c__[1];
        }
        i__1 = *ns - 1;
        for (k = 1; k <= i__1; ++k) {
          fact = 1. / (weight_1.c__[*ns - k] - weight_1.c__[*ns - k + 1]);
          i__2 = *n;
          for (i__ = 1; i__ <= i__2; ++i__) {
            cont[i__ + k * *n] =
                (zz[i__ + (*ns - k - 1) * *n] - zz[i__ + (*ns - k) * *n]) *
                fact;
          }
        }
        i__1 = *ns;
        for (j = 2; j <= i__1; ++j) {
          i__2 = j;
          for (k = *ns; k >= i__2; --k) {
            fact = 1. / (weight_1.c__[*ns - k] - weight_1.c__[*ns - k + j]);
            i__3 = *n;
            for (i__ = 1; i__ <= i__3; ++i__) {
              cont[i__ + k * *n] =
                  (cont[i__ + k * *n] - cont[i__ + (k - 1) * *n]) * fact;
            }
          }
        }
      }
      /* *** *** *** *** *** *** *** */
      /* *** *** *** *** *** *** *** */
    } else {
      if (*ns == 7) {
        /* *** *** *** *** *** *** *** */
        /* *** *** *** *** *** *** *** */
        if (first || *startn || change) {
          i__1 = nns;
          for (i__ = 1; i__ <= i__1; ++i__) {
            zz[i__] = 0.;
            ff[i__] = 0.;
          }
        } else {
          hquot = *h__ / hold;
          i__1 = *ns;
          for (k = 1; k <= i__1; ++k) {
            ccq = weight_1.c__[k] * hquot;
            i__2 = *n;
            for (i__ = 1; i__ <= i__2; ++i__) {
              val = cont[i__ + *ns * *n];
              for (l = *ns - 1; l >= 1; --l) {
                val = cont[i__ + l * *n] +
                      (ccq - weight_1.c__[*ns - l] + 1.) * val;
              }
              zz[i__ + (k - 1) * *n] = ccq * val;
            }
          }
          i__1 = *n;
          for (i__ = 1; i__ <= i__1; ++i__) {
            z1i = zz[i__];
            z2i = zz[i__ + *n];
            z3i = zz[i__ + n2];
            z4i = zz[i__ + n3];
            z5i = zz[i__ + n4];
            z6i = zz[i__ + n5];
            z7i = zz[i__ + n6];
            ff[i__] = coe7.ti711 * z1i + coe7.ti712 * z2i + coe7.ti713 * z3i +
                      coe7.ti714 * z4i + coe7.ti715 * z5i + coe7.ti716 * z6i +
                      coe7.ti717 * z7i;
            ff[i__ + *n] = coe7.ti721 * z1i + coe7.ti722 * z2i +
                           coe7.ti723 * z3i + coe7.ti724 * z4i +
                           coe7.ti725 * z5i + coe7.ti726 * z6i +
                           coe7.ti727 * z7i;
            ff[i__ + n2] = coe7.ti731 * z1i + coe7.ti732 * z2i +
                           coe7.ti733 * z3i + coe7.ti734 * z4i +
                           coe7.ti735 * z5i + coe7.ti736 * z6i +
                           coe7.ti737 * z7i;
            ff[i__ + n3] = coe7.ti741 * z1i + coe7.ti742 * z2i +
                           coe7.ti743 * z3i + coe7.ti744 * z4i +
                           coe7.ti745 * z5i + coe7.ti746 * z6i +
                           coe7.ti747 * z7i;
            ff[i__ + n4] = coe7.ti751 * z1i + coe7.ti752 * z2i +
                           coe7.ti753 * z3i + coe7.ti754 * z4i +
                           coe7.ti755 * z5i + coe7.ti756 * z6i +
                           coe7.ti757 * z7i;
            ff[i__ + n5] = coe7.ti761 * z1i + coe7.ti762 * z2i +
                           coe7.ti763 * z3i + coe7.ti764 * z4i +
                           coe7.ti765 * z5i + coe7.ti766 * z6i +
                           coe7.ti767 * z7i;
            ff[i__ + n6] = coe7.ti771 * z1i + coe7.ti772 * z2i +
                           coe7.ti773 * z3i + coe7.ti774 * z4i +
                           coe7.ti775 * z5i + coe7.ti776 * z6i +
                           coe7.ti777 * z7i;
          }
        }
        /* *** *** *** *** *** *** *** */
        /*  LOOP FOR THE SIMPLIFIED NEWTON ITERATION */
        /* *** *** *** *** *** *** *** */
        newt = 0;
        nit = *nit1 + 10;
        expmi = 1. / expmns;
        /* Computing std::max */
        /* Computing std::min */
        d__5 = expmi - 1.;
        d__3 = .03, d__4 = pow(rtol1, d__5);
        d__1 = *uround * 10 / rtol1, d__2 = std::min(d__3, d__4);
        fnewt = std::max(d__1, d__2);
        d__1 = std::max(faccon, *uround);
        faccon = pow(d__1, c_b100);
        theta = std::abs(*thet);
      L240:
        if (newt >= nit) {
          goto L78;
        }
        /* ---     COMPUTE THE RIGHT-HAND SIDE */
        i__1 = *ns - 1;
        for (k = 0; k <= i__1; ++k) {
          iadd = k * *n;
          i__2 = *n;
          for (i__ = 1; i__ <= i__2; ++i__) {
            cont[i__] = y[i__] + zz[iadd + i__];
          }
          d__1 = *x + weight_1.c__[k + 1] * *h__;
          fcn(n, &d__1, &cont[1], &zz[iadd + 1]);
        }
        *nfcn += *ns;
        /* ---     SOLVE THE LINEAR SYSTEMS */
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
          z1i = zz[i__];
          z2i = zz[i__ + *n];
          z3i = zz[i__ + n2];
          z4i = zz[i__ + n3];
          z5i = zz[i__ + n4];
          z6i = zz[i__ + n5];
          z7i = zz[i__ + n6];
          zz[i__] = coe7.ti711 * z1i + coe7.ti712 * z2i + coe7.ti713 * z3i +
                    coe7.ti714 * z4i + coe7.ti715 * z5i + coe7.ti716 * z6i +
                    coe7.ti717 * z7i;
          zz[i__ + *n] = coe7.ti721 * z1i + coe7.ti722 * z2i +
                         coe7.ti723 * z3i + coe7.ti724 * z4i +
                         coe7.ti725 * z5i + coe7.ti726 * z6i + coe7.ti727 * z7i;
          zz[i__ + n2] = coe7.ti731 * z1i + coe7.ti732 * z2i +
                         coe7.ti733 * z3i + coe7.ti734 * z4i +
                         coe7.ti735 * z5i + coe7.ti736 * z6i + coe7.ti737 * z7i;
          zz[i__ + n3] = coe7.ti741 * z1i + coe7.ti742 * z2i +
                         coe7.ti743 * z3i + coe7.ti744 * z4i +
                         coe7.ti745 * z5i + coe7.ti746 * z6i + coe7.ti747 * z7i;
          zz[i__ + n4] = coe7.ti751 * z1i + coe7.ti752 * z2i +
                         coe7.ti753 * z3i + coe7.ti754 * z4i +
                         coe7.ti755 * z5i + coe7.ti756 * z6i + coe7.ti757 * z7i;
          zz[i__ + n5] = coe7.ti761 * z1i + coe7.ti762 * z2i +
                         coe7.ti763 * z3i + coe7.ti764 * z4i +
                         coe7.ti765 * z5i + coe7.ti766 * z6i + coe7.ti767 * z7i;
          zz[i__ + n6] = coe7.ti771 * z1i + coe7.ti772 * z2i +
                         coe7.ti773 * z3i + coe7.ti774 * z4i +
                         coe7.ti775 * z5i + coe7.ti776 * z6i + coe7.ti777 * z7i;
        }
        slvrar(n, &fjac[fjac_offset], ldjac, mljac, mujac, &fmas[fmas_offset],
               ldmas, mlmas, mumas, m1, m2, nm1, &fac1, &e1[e1_offset], lde1,
               &zz[1], &ff[1], &ip1[1], &iphes[1], &ier, ijob, linal_1);
        for (k = 1; k <= 3; ++k) {
          // v--- this was: iad = (k - 1 << 1) * *nm1 + 1;
          iad = ((k - 1) << 1) * *nm1 + 1;
          slvrai(n, &fjac[fjac_offset], ldjac, mljac, mujac, &fmas[fmas_offset],
                 ldmas, mlmas, mumas, m1, m2, nm1, &alphn[k - 1], &betan[k - 1],
                 &ee2[iad * ee2_dim1 + 1], &ee2[(iad + *nm1) * ee2_dim1 + 1],
                 lde1, &zz[((k << 1) - 1) * *n + 1], &zz[(k << 1) * *n + 1],
                 &ff[((k << 1) - 1) * *n + 1], &ff[(k << 1) * *n + 1], &cont[1],
                 &ip2[(k - 1) * *nm1 + 1], &iphes[1], &ier, ijob, linal_1);
        }
        ++(*nsol);
        ++newt;
        dyno = 0.;
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
          denom = scal[i__];
          i__2 = *ns - 1;
          for (k = 0; k <= i__2; ++k) {
            /* Computing 2nd power */
            d__1 = zz[i__ + k * *n] / denom;
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
            dyth = faccon * dyno * pow(theta, i__1) / fnewt;
            if (dyth >= 1.) {
              /* Computing std::max */
              d__1 = 1e-4, d__2 = std::min(20., dyth);
              qnewt = std::max(d__1, d__2);
              d__1 = -1. / (nit + 4. - 1 - newt);
              hhfac = pow(qnewt, d__1) * 0.8;
              *h__ = hhfac * *h__;
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
        dynold = std::max(dyno, *uround);
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
          z1i = ff[i__] + zz[i__];
          z2i = ff[i__ + *n] + zz[i__ + *n];
          z3i = ff[i__ + n2] + zz[i__ + n2];
          z4i = ff[i__ + n3] + zz[i__ + n3];
          z5i = ff[i__ + n4] + zz[i__ + n4];
          z6i = ff[i__ + n5] + zz[i__ + n5];
          z7i = ff[i__ + n6] + zz[i__ + n6];
          ff[i__] = z1i;
          ff[i__ + *n] = z2i;
          ff[i__ + n2] = z3i;
          ff[i__ + n3] = z4i;
          ff[i__ + n4] = z5i;
          ff[i__ + n5] = z6i;
          ff[i__ + n6] = z7i;
          zz[i__] = coe7.t711 * z1i + coe7.t712 * z2i + coe7.t713 * z3i +
                    coe7.t714 * z4i + coe7.t715 * z5i + coe7.t716 * z6i +
                    coe7.t717 * z7i;
          zz[i__ + *n] = coe7.t721 * z1i + coe7.t722 * z2i + coe7.t723 * z3i +
                         coe7.t724 * z4i + coe7.t725 * z5i + coe7.t726 * z6i +
                         coe7.t727 * z7i;
          zz[i__ + n2] = coe7.t731 * z1i + coe7.t732 * z2i + coe7.t733 * z3i +
                         coe7.t734 * z4i + coe7.t735 * z5i + coe7.t736 * z6i +
                         coe7.t737 * z7i;
          zz[i__ + n3] = coe7.t741 * z1i + coe7.t742 * z2i + coe7.t743 * z3i +
                         coe7.t744 * z4i + coe7.t745 * z5i + coe7.t746 * z6i +
                         coe7.t747 * z7i;
          zz[i__ + n4] = coe7.t751 * z1i + coe7.t752 * z2i + coe7.t753 * z3i +
                         coe7.t754 * z4i + coe7.t755 * z5i + coe7.t756 * z6i +
                         coe7.t757 * z7i;
          zz[i__ + n5] = coe7.t761 * z1i + coe7.t762 * z2i + coe7.t763 * z3i +
                         coe7.t764 * z4i + coe7.t765 * z5i + coe7.t766 * z6i +
                         coe7.t767 * z7i;
          zz[i__ + n6] = coe7.t771 * z1i + z2i + z4i + z6i;
        }
        if (faccon * dyno > fnewt) {
          goto L240;
        }
        /* --- ERROR ESTIMATION */
        estrav(n, &fjac[fjac_offset], ldjac, mljac, mujac, &fmas[fmas_offset],
               ldmas, mlmas, mumas, h__, dd, fcn, nfcn, &y0[1], &y[1], ijob, x,
               m1, m2, nm1, ns, &nns, &e1[e1_offset], lde1, &zz[1], &cont[1],
               &ff[1], &ip1[1], &iphes[1], &scal[1], &err, &first, &reject,
               &fac1, linal_1);
        /*       --- COMPUTE FINITE DIFFERENCES FOR DENSE OUTPUT */
        if (err < 1.) {
          i__1 = *n;
          for (i__ = 1; i__ <= i__1; ++i__) {
            y[i__] += zz[i__ + (*ns - 1) * *n];
            cont[i__ + *ns * *n] = zz[i__] / weight_1.c__[1];
          }
          i__1 = *ns - 1;
          for (k = 1; k <= i__1; ++k) {
            fact = 1. / (weight_1.c__[*ns - k] - weight_1.c__[*ns - k + 1]);
            i__2 = *n;
            for (i__ = 1; i__ <= i__2; ++i__) {
              cont[i__ + k * *n] =
                  (zz[i__ + (*ns - k - 1) * *n] - zz[i__ + (*ns - k) * *n]) *
                  fact;
            }
          }
          i__1 = *ns;
          for (j = 2; j <= i__1; ++j) {
            i__2 = j;
            for (k = *ns; k >= i__2; --k) {
              fact = 1. / (weight_1.c__[*ns - k] - weight_1.c__[*ns - k + j]);
              i__3 = *n;
              for (i__ = 1; i__ <= i__3; ++i__) {
                cont[i__ + k * *n] =
                    (cont[i__ + k * *n] - cont[i__ + (k - 1) * *n]) * fact;
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
          i__1 = *ns;
          for (i__ = 1; i__ <= i__1; ++i__) {
            zz[i__] = 0.;
            ff[i__] = 0.;
          }
        } else {
          hquot = *h__ / hold;
          i__1 = *n;
          for (i__ = 1; i__ <= i__1; ++i__) {
            z1i = hquot * cont[i__ + *n];
            zz[i__] = z1i;
            ff[i__] = z1i;
          }
        }
        /* *** *** *** *** *** *** *** */
        /*  LOOP FOR THE SIMPLIFIED NEWTON ITERATION */
        /* *** *** *** *** *** *** *** */
        newt = 0;
        nit = *nit1 - 3;
        expmi = 1. / expmns;
        /* Computing std::max */
        d__1 = *uround * 10 / rtol1;
        fnewt = std::max(d__1, .03);
        d__1 = std::max(faccon, *uround);
        faccon = pow(d__1, c_b100);
        theta = std::abs(*thet);
      L440:
        if (newt >= nit) {
          goto L78;
        }
        /* ---     COMPUTE THE RIGHT-HAND SIDE */
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
          cont[i__] = y[i__] + zz[i__];
        }
        fcn(n, &xph, &cont[1], &zz[1]);
        ++(*nfcn);
        /* ---     SOLVE THE LINEAR SYSTEMS */
        slvrar(n, &fjac[fjac_offset], ldjac, mljac, mujac, &fmas[fmas_offset],
               ldmas, mlmas, mumas, m1, m2, nm1, &fac1, &e1[e1_offset], lde1,
               &zz[1], &ff[1], &ip1[1], &iphes[1], &ier, ijob, linal_1);
        ++(*nsol);
        ++newt;
        dyno = 0.;
        i__1 = *n;
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
            dyth = faccon * dyno * pow(theta, i__1) / fnewt;
            if (dyth >= 1.) {
              /* Computing std::max */
              d__1 = 1e-4, d__2 = std::min(20., dyth);
              qnewt = std::max(d__1, d__2);
              d__1 = -1. / (nit + 4. - 1 - newt);
              hhfac = pow(qnewt, d__1) * .8;
              *h__ = hhfac * *h__;
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
        dynold = std::max(dyno, *uround);
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
          f1i = ff[i__] + zz[i__];
          ff[i__] = f1i;
          zz[i__] = f1i;
        }
        if (faccon * dyno > fnewt) {
          goto L440;
        }
        /* --- ERROR ESTIMATION */
        estrav(n, &fjac[fjac_offset], ldjac, mljac, mujac, &fmas[fmas_offset],
               ldmas, mlmas, mumas, h__, dd, fcn, nfcn, &y0[1], &y[1], ijob, x,
               m1, m2, nm1, ns, &nns, &e1[e1_offset], lde1, &zz[1], &cont[1],
               &ff[1], &ip1[1], &iphes[1], &scal[1], &err, &first, &reject,
               &fac1, linal_1);
        /*       --- COMPUTE FINITE DIFFERENCES FOR DENSE OUTPUT */
        if (err < 1.) {
          i__1 = *n;
          for (i__ = 1; i__ <= i__1; ++i__) {
            y[i__] += zz[i__];
            cont[i__ + *n] = zz[i__];
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
  /* Computing std::min */
  d__1 = *safe, d__2 = ((nit << 1) + 1) * *safe / (newt + (nit << 1));
  fac = std::min(d__1, d__2);
  /* Computing std::max */
  /* Computing std::min */
  d__3 = *facl, d__4 = pow(err, expo) / fac;
  d__1 = *facr, d__2 = std::min(d__3, d__4);
  quot = std::max(d__1, d__2);
  hnew = *h__ / quot;
  /* *** *** *** *** *** *** *** */
  /*  IS THE ERROR SMALL ENOUGH ? */
  /* *** *** *** *** *** *** *** */
  if (err < 1.) {
    /* --- STEP IS ACCEPTED */
    first = false;
    ++(*naccpt);
    if (*pred && !change) {
      /*       --- PREDICTIVE CONTROLLER OF GUSTAFSSON */
      if (*naccpt > 1) {
        /* Computing 2nd power */
        d__2 = err;
        d__1 = d__2 * d__2 / erracc;
        facgus = hacc / *h__ * pow(d__1, expo) / *safe;
        /* Computing std::max */
        d__1 = *facr, d__2 = std::min(*facl, facgus);
        facgus = std::max(d__1, d__2);
        quot = std::max(quot, facgus);
        hnew = *h__ / quot;
      }
      hacc = *h__;
      erracc = std::max(.01, err);
    }
    xold = *x;
    hold = *h__;
    *x = xph;
    /*       --- UPDATE SCALING */
    if (*itol == 0) {
      i__1 = *n;
      for (i__ = 1; i__ <= i__1; ++i__) {
        scal[i__] = atol1 + rtol1 * (d__1 = y[i__], std::abs(d__1));
      }
    } else {
      i__1 = *n;
      for (i__ = 1; i__ <= i__1; ++i__) {
        quott = atol[i__] / rtol[i__];
        rtol1 = pow(rtol[i__], expmns) * .1;
        atol1 = rtol1 * quott;
        scal[i__] = atol1 + rtol1 * (d__1 = y[i__], std::abs(d__1));
      }
    }
    if (*iout != 0) {
      nrsol = *naccpt + 1;
      weight_1.xsol = *x;
      xosol = xold;
      i__1 = *n;
      for (i__ = 1; i__ <= i__1; ++i__) {
        cont[i__] = y[i__];
      }
      nsolu = *n;
      weight_1.hsol = hold;
      solout(&nrsol, &xosol, &weight_1.xsol, &y[1], &cont[1], &lrc, &nsolu,
             &irtrn, weight_1);
      if (irtrn < 0) {
        goto L179;
      }
    }
    caljac = false;
    if (last) {
      *h__ = hopt;
      *idid = 1;
      return 0;
    }
    fcn(n, x, &y[1], &y0[1]);
    ++(*nfcn);
    /* Computing std::min */
    d__1 = std::abs(hnew);
    hnew = posneg * std::min(d__1, hmaxn);
    hopt = hnew;
    hopt = std::min(*h__, hnew);
    if (reject) {
      /* Computing std::min */
      d__1 = std::abs(hnew), d__2 = std::abs(*h__);
      hnew = posneg * std::min(d__1, d__2);
    }
    reject = false;
    if ((*x + hnew / *quot1 - *xend) * posneg >= 0.) {
      *h__ = *xend - *x;
      last = true;
    } else {
      qt = hnew / *h__;
      hhfac = *h__;
      if (theta <= *thet && qt >= *quot1 && qt <= *quot2) {
        ikeep = 1;
        goto L30;
      }
      *h__ = hnew;
    }
    hhfac = *h__;
    if (theta <= *thet) {
      goto L20;
    }
    goto L10;
  } else {
    /* --- STEP IS REJECTED */
    reject = true;
    last = false;
    if (first) {
      *h__ *= .1;
      hhfac = .1;
    } else {
      hhfac = hnew / *h__;
      *h__ = hnew;
    }
    if (*naccpt >= 1) {
      ++(*nrejct);
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
  *h__ *= .5;
  hhfac = .5;
  reject = true;
  last = false;
  if (caljac) {
    goto L20;
  }
  goto L10;
  /* --- FAIL EXIT */
L176:
  std::cerr << "Encountered repeatedly singular matrix, ier = " << ier << "\n";
  *idid = -4;
  return 0;
L177:
  std::cerr << "Step-size too small: h = " << *h__ << "\n";
  *idid = -3;
  return 0;
L178:
  std::cerr << "More than " << *nmax << " steps needed\n";
  *idid = -2;
  return 0;
  /* --- EXIT CAUSED BY SOLOUT */
L179:
  *idid = 2;
  return 0;
}

//===========================================================================
//---- Radau Integrator Interface -------------------------------------------
//===========================================================================

int radau(int *n, F_fcn fcn, double *x, double *y, double *xend, double *h__,
          double *rtol, double *atol, int *itol, F_jac jac, int *ijac,
          int *mljac, int *mujac, F_mas mas, int *imas, int *mlmas, int *mumas,
          F_solout solout, int *iout, double *work, int *lwork, int *iwork,
          int *liwork, int *idid) {
  /* System generated locals */
  int i__1;

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

  /* Function Body */
  nfcn = 0;
  njac = 0;
  nstep = 0;
  naccpt = 0;
  nrejct = 0;
  ndec = 0;
  nsol = 0;
  arret = false;
  /* -------- NUMBER std::maxIMAL AND std::minIMAL OF STAGES  NS */
  if (iwork[11] == 0) {
    nsmin = 3;
  } else {
    nsmin = std::max(1, iwork[11]);
    if (iwork[11] >= 2) {
      nsmin = std::max(3, iwork[11]);
    }
    if (iwork[11] >= 4) {
      nsmin = std::max(5, iwork[11]);
    }
    if (iwork[11] >= 6) {
      nsmin = 7;
    }
  }
  if (iwork[12] == 0) {
    nsmax = 7;
  } else {
    nsmax = std::min(7, iwork[12]);
    if (iwork[12] <= 6) {
      nsmax = std::min(5, iwork[12]);
    }
    if (iwork[12] <= 4) {
      nsmax = std::min(3, iwork[12]);
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
      std::cerr << "Invalid input: iwork[12] (aka nsus) = " << iwork[13]
                << "\n";
      arret = true;
    }
  }
  /* -------- nmax , THE std::maxIMAL NUMBER OF STEPS ----- */
  if (iwork[2] == 0) {
    nmax = 100000;
  } else {
    nmax = iwork[2];
    if (nmax <= 0) {
      std::cerr << "Invalid input: iwork[1] (aka nmax) = " << iwork[2] << "\n";
      arret = true;
    }
  }
  /* -------- NIT    std::maxIMAL NUMBER OF NEWTON ITERATIONS */
  if (iwork[3] == 0) {
    nit = 7;
  } else {
    nit = iwork[3];
    if (nit <= 0 || nit > 50) {
      std::cerr << "Invalid input: iwork[2] (aka nit) = " << iwork[3] << "\n";
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
    nind1 = *n;
  }
  if (nind1 + nind2 + nind3 != *n) {
    std::cerr
        << "Invalid input: iwork[4],iwork[5] or iwork[6] (aka nind1,2,3) = "
        << nind1 << ", " << nind2 << ", " << nind3 << "\n";
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
  nm1 = *n - m1;
  if (m1 == 0) {
    m2 = *n;
  }
  if (m2 == 0) {
    m2 = m1;
  }
  if (m1 < 0 || m2 < 0 || m1 + m2 > *n) {
    std::cerr << "Invalid input: iwork[8] or iwork[9] (aka m1,m2) = " << m1
              << ", " << m2 << "\n";
    arret = true;
  }
  /* -------- UROUND   SMALLEST NUMBER SATISFYING 1.0D0+UROUND>1.0D0 */
  if (work[1] == 0.) {
    uround = 1e-16;
  } else {
    uround = work[1];
    if (uround <= 1e-19 || uround >= 1.) {
      std::cerr << "Coefficients have 20 digits, uround = " << work[1] << "\n";
      arret = true;
    }
  }
  /* --------- CHECK IF TOLERANCES ARE O.K. */
  if (*itol == 0) {
    if (atol[1] <= 0. || rtol[1] <= uround * 10.) {
      std::cerr << "Tolerances are too small\n";
      arret = true;
    }
  } else {
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      if (atol[i__] <= 0. || rtol[i__] <= uround * 10.) {
        std::cerr << "Tolerances " << i__ << " are too small\n";
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
      std::cerr << "Invalid input: work[1] (aka safe) = " << work[2] << "\n";
      arret = true;
    }
  }
  /* ------ THET     DECIDES WHETHER THE JACOBIAN SHOULD BE RECOMPUTED; */
  if (work[3] == 0.) {
    thet = .001;
  } else {
    thet = work[3];
    if (thet >= 1.) {
      std::cerr << "Invalid input: work[2] (aka thet) = " << work[3] << "\n";
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
    std::cerr << "Invalid input: work[4] or work[5] (aka quot1,2) = " << quot1
              << ", " << quot2 << "\n";
    arret = true;
  }
  /* -------- std::maxIMAL STEP SIZE */
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
    std::cerr << "Invalid input: work[7] or work[7] (aka facl,r) = " << work[9]
              << ", " << work[9] << "\n";
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
      std::cerr
          << "Bandwidth of 'mas' is not smaller than bandwidth of 'jac'\n";
      arret = true;
    }
  } else {
    ldmas = 0;
    if (jband) {
      ijob = 2;
    } else {
      ijob = 1;
      if (*n > 2 && iwork[1] != 0) {
        ijob = 7;
      }
    }
  }
  ldmas2 = std::max(1, ldmas);
  /* ------ HESSENBERG OPTION ONLY FOR EXPLICIT EQU. WITH FULL JACOBIAN */
  if ((implct || jband) && ijob == 7) {
    std::cerr
        << "Hessenberg option only for explicit equations with full jacobian\n";
    arret = true;
  }
  /* ------- PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN WORK ----- */
  nns = ns * *n;
  nm1ns = ns * nm1;
  nmee = (ns - 1) * nm1;
  iezz = 21;
  iey0 = iezz + nns;
  iescal = iey0 + *n;
  ieff = iescal + *n;
  iecon = ieff + nns;
  iejac = iecon + nns + *n;
  iemas = iejac + *n * ldjac;
  iee1 = iemas + nm1 * ldmas;
  iee = iee1 + nm1 * lde1;
  /* ------ TOTAL STORAGE REQUIREMENT ----------- */
  istore = iee + nmee * lde1 - 1;
  if (istore > *lwork) {
    std::cerr << "Array 'work' is too small = " << istore << "\n";
    arret = true;
  }
  /* ------- ENTRY POINTS FOR int WORKSPACE ----- */
  ieip1 = 21;
  ieip2 = ieip1 + nm1;
  ieiph = ieip2 + nm1 * (ns - 1) / 2;
  /* --------- TOTAL REQUIREMENT --------------- */
  istore = ieiph + nm1 - 1;
  if (istore > *liwork) {
    std::cerr << "Array 'iwork' is too small = " << istore << "\n";
    arret = true;
  }
  /* ------ WHEN A FAIL HAS OCCURED, WE RETURN WITH IDID=-1 */
  if (arret) {
    *idid = -1;
    return 0;
  }
  /* -------- CALL TO CORE INTEGRATOR ------------ */
  radcov(n, fcn, x, &y[1], xend, &hmax, h__, &rtol[1], &atol[1], itol, &nsus,
         jac, ijac, mljac, mujac, mas, mlmas, mumas, solout, iout, idid, &nmax,
         &uround, &safe, &thet, &quot1, &quot2, &nit, &ijob, &startn, &nind1,
         &nind2, &nind3, &pred, &facl, &facr, &m1, &m2, &nm1, &nsmin, &ns, &nns,
         &nm1ns, &nmee, &implct, &jband, &ldjac, &lde1, &ldmas2, &work[iezz],
         &work[iey0], &work[iescal], &work[ieff], &work[iejac], &work[iee1],
         &work[iee], &work[iemas], &work[iecon], &iwork[ieip1], &iwork[ieip2],
         &iwork[ieiph], &vitu, &vitd, &hhou, &hhod, &nfcn, &njac, &nstep,
         &naccpt, &nrejct, &ndec, &nsol);
  iwork[13] = nsus;
  iwork[14] = nfcn;
  iwork[15] = njac;
  iwork[16] = nstep;
  iwork[17] = naccpt;
  iwork[18] = nrejct;
  iwork[19] = ndec;
  iwork[20] = nsol;
  return 0;
}

} // namespace stiff

#endif // STIFF_RADAU_HPP
