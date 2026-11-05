/**
 * my_diag.c - diagnostics of the 3-shunt ADC sampling point (MCSDK R3_2_G4xx).
 *
 * READ-ONLY: nothing in the MCSDK library is modified. The module only looks at
 * the PWM handle and TIM1->CCR4 and counts which sampling case is in use:
 *
 *   case 1 : CCR4 == Half_PWMPeriod-1, sampling in the middle of the PWM period.
 *            TRISE/TNOISE/dead time have NO effect here.
 *   case 2/3: CCR4 moved (or falling edge trigger used) - window of low side
 *            is shorter than Tafter. This is where TRISE/TNOISE act.
 *   est    : useEstCurrent == true (the current is not measured, it is estimated)
 *
 * USE
 *   1. Add my_diag.c to the project (Src/).
 *   2. Call MyDiag_Init() once after MC_boot (e.g. where myapp init runs).
 *   3. Call MyDiag_Step() as often as possible. Minimum: from the Medium
 *      Frequency Task or myapp periodic function (1 kHz gives a statistic,
 *      not every PWM cycle). Best: from a TIM1 update/ADC JEOC hook if you
 *      already have one (then every PWM cycle is seen).
 *   4. Watch g_myDiag in the debugger / Pilot (all fields volatile, 32-bit) or
 *      call MyDiag_Report(buf, size) once per second to print a line.
 *   5. MyDiag_Reset() clears counters (call it before each test run, e.g. at
 *      motor start; min/max are cleared too).
 *
 * Fields in g_myDiag:
 *   n_total, n_case1, n_case23, n_est       - counters
 *   pct_case23_x10, pct_est_x10             - percentage * 10 over the whole run
 *   headroom_min                            - min of (Half_PWMPeriod-lowDuty-Tafter), timer ticks
 *                                             (<=0 means case 1 is lost)
 *   headroom_ns_min                         - same, nanoseconds
 *   lowDuty_max, duty_pct_x10_max           - max lowDuty and its % of Half_PWMPeriod, x10
 *   vs_s16_max                              - max |Vqd| in s16 digits (full scale 32767)
 *   vs_pct_x10_max                          - same, % of 32767 x10
 *   iq_abs_max_s16                          - max |Iq| (s16 digit) when case23 was seen
 *   iq_at_case23                            - last Iq (s16 digit) when case23 was seen
 *   vs_at_case23                            - last |Vqd| when case23 was seen
 *   speed_at_case23_rpm                     - last speed (0.1 Hz units->rpm) when case23 was seen
 *   Tafter, Half, Tbefore, Tcase2, Tcase3, Tsampling - the values that are really in the handle
 *
 * v12: PASSIVE IDENTIFICATION (see the block "PASSIVE IDENTIFICATION" below). No new calls are needed:
 *   MyDiag_Step() collects windows, MyDiag_CycleEnd() starts the computation, the result appears in
 *   g_myDiagIdLast (and g_myDiagIdHist[same index as g_myDiagHist]) a few ms later, in the Stop state.
 *   Live Expressions: g_myDiagIdLast.fit[0..3], .L_bin_uH_x10[], .n_bin[], .shift_now_x10, g_myDiagIdAll[0..3].
 *   Also new in g_myDiagHist: idref_mean_s16 (mean Id reference in the same window as id_mean_s16).
 *
 * TEST PLAN FOR TRISE
 *   run the same cycle with TRISE_NS=200, then 2000, then 5000. Compare
 *   pct_case23_x10, headroom_min, vs_pct_x10_max. If pct_case23 stays ~0 in all runs
 *   the parameter simply never gets used in your operating range.
 */

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "main.h"
#include "mc_type.h"
#include "mc_config.h"
#include "mc_interface.h"
#include "pqd_motor_power_measurement.h"
#include "pwm_curr_fdbk.h"
#include "r3_2_g4xx_pwm_curr_fdbk.h"
#include "pmsm_motor_parameters.h"   /* POLE_PAIR_NUM */
#include "hall_speed_pos_fdbk.h"     /* HALL_M1.PhaseShift (read only) */

/* ---------- configuration ---------- */
#define MYDIAG_MOTOR_IDX     0u           /* M1 */
#define MYDIAG_TIM_CLK_MHZ   170u         /* ADV_TIM_CLK_MHz of the project */
#define MYDIAG_NEG_RPM       100          /* speed below -100 rpm counts as reverse rotation */
#define MYDIAG_ANG_RPM_MIN   2500         /* angle estimate only above this speed (EMF dominates) */
#define MYDIAG_ANG_IQ_MAX    1400         /* ... and below this |Iq| in s16 digits (about 15 A, scale not verified) */
#define MYDIAG_ANG_STEADY_RPM 150         /* window only when |speed - speed_ref| < this (steady state) */
#define MYDIAG_MIN_SAMPLES   200u         /* shorter runs are ignored by MyDiag_CycleEnd */

/* --- passive identification of R, L, lambda(Ke) and the frame angle error (see the block below) --- */
#define MYDIAG_ID_ENABLE      1
#define MYDIAG_POLE_PAIRS     POLE_PAIR_NUM   /* 4 BLACK, 7 DOMEL (pmsm_motor_parameters.h) */
#define MYDIAG_DIGITS_PER_A   91.5f           /* s16 digits of Id/Iq per 1 A (1 mOhm shunt, gain 4.606, 3.3 V, 16 bit) */
#define MYDIAG_VQD_FS_DIV     1.7320508f      /* Vd/Vq full scale 32767 = Vbus / this (SVPWM linear limit) */
#define MYDIAG_ID_WIN_MS      64u             /* window length = number of Step() calls (Step must run at 1 kHz) */
#define MYDIAG_ID_NWIN        192u            /* ring of windows: 12 bytes each (192 -> 2.3 kB, about 12 s of data) */
#define MYDIAG_ID_RPM_MIN     800             /* ignore slower windows (Hall interpolation, no EMF) */
#define MYDIAG_ID_SPAN_PCT    8               /* window rejected if max-min speed inside it > this % of the mean */
#define MYDIAG_ID_STEADY_RPM  300             /* window rejected if |speed - speed_ref| > this at any ms */
#define MYDIAG_ID_VS_MAX_X10  900             /* window rejected if |Vqd| > 90.0 % of full scale (voltage clipped) */
#define MYDIAG_ID_MIN_WIN     12u             /* fewer windows -> fit not valid */
#define MYDIAG_ID_TAU_SPAN_RPM 600            /* 'tau' fit only if the speed range of the windows is wider than this */
#define MYDIAG_ID_REF_RPM     3500            /* speed for which delta_ref is extrapolated (working speed) */
#define MYDIAG_USE_HALL_M1    1               /* read the current HALL_M1.PhaseShift for the report */
#define MYDIAG_USE_IDREF      1               /* log mean Id reference (MCI_GetIqdref) */

/* --- Hall sector diagnostics (widths of the 6 sectors, predicted angle ripple vs HALL_AVERAGING_FIFO_DEPTH) --- */
#define MYDIAG_HALL_ENABLE    1
#define MYDIAG_HALL_RPM_MIN   1500            /* collect only above this speed ... */
#define MYDIAG_HALL_STEADY_RPM 150            /* ... and when |speed - speed_ref| < this */
#define MYDIAG_HALL_REV_TOL_PCT 3             /* a revolution is used only if its period differs from the previous one by < this % */
/* ----------------------------------- */

typedef struct
{
  volatile uint32_t n_total;
  volatile uint32_t n_case1;
  volatile uint32_t n_case23;
  volatile uint32_t n_est;
  volatile uint32_t pct_case23_x10;
  volatile uint32_t pct_est_x10;

  volatile int32_t  headroom_min;
  volatile int32_t  headroom_ns_min;
  volatile uint32_t lowDuty_max;
  volatile uint32_t duty_pct_x10_max;

  volatile uint32_t vs_s16_max;
  volatile uint32_t vs_pct_x10_max;
  volatile uint32_t imag_max_c1;      /* max |I| while sampled in the middle (case 1) */
  volatile uint32_t imag_max_c23;     /* max |I| while sample point was moved (case 2/3) */
  volatile uint32_t imag_t_ms;        /* ms since run start at the |I| maximum */
  volatile int32_t  imag_rpm;         /* speed at the |I| maximum */
  volatile uint32_t i_mag_max_s16;    /* max sqrt(Id^2+Iq^2), digits */
  volatile uint32_t id_abs_max_s16;   /* max |Id|, digits */
  volatile uint32_t iq_abs_max_s16;
  volatile int32_t  iq_at_case23;
  volatile uint32_t vs_at_case23;
  volatile int32_t  speed_at_case23_rpm;

  volatile int32_t  speed_min_rpm;     /* signed, negative = shaft turns backwards */
  volatile int32_t  speed_max_rpm;
  volatile uint32_t n_neg;             /* samples (ms) with speed < -MYDIAG_NEG_RPM */

  /* Frame misalignment estimate. At light load and constant speed the voltage vector is
   * the back-EMF: with a correct rotor angle it lies on q (Vd ~ 0). The angle of the mean
   * vector = atan2(Vd, Vq) is the angle error (sign shows the direction to move HALL_PHASE_SHIFT). */
  volatile int32_t  sum_vd;
  volatile int32_t  sum_vq;
  volatile uint32_t n_ang;
  volatile float    sum_iq2;           /* sum of squares, for the ripple (std) of Iq / Id */
  volatile float    sum_id2;
  volatile float    sum_rpm;           /* for the mean speed / bus voltage of the window (Ls estimate) */
  volatile float    sum_vbus;
  volatile int32_t  sum_iq;
  volatile int32_t  sum_id;
  volatile int32_t  sum_idref;         /* sum of the Id reference, same window */
  volatile float    sum_pw;            /* sum of PQD average electrical power [W], same window as sum_vd */

  volatile uint32_t Tafter;
  volatile uint32_t Half;
  volatile uint32_t Tbefore;
  volatile uint32_t Tcase2;
  volatile uint32_t Tcase3;
  volatile uint32_t Tsampling;
} MyDiag_t;

volatile MyDiag_t g_myDiag;

/* Summary of the last cycles (ring of 8). Filled by MyDiag_CycleEnd(). */
typedef struct
{
  uint32_t n_total;
  uint16_t pct_case23_x10;   /* % of PWM-cycle samples not in the middle, x10 */
  uint16_t pct_est_x10;      /* % with estimated current, x10 */
  int16_t  headroom_ns_min;  /* min window margin vs Tafter, ns (<=0: case 1 lost) */
  uint16_t duty_pct_x10_max; /* max lowDuty / Half_PWMPeriod, % x10 */
  uint16_t vs_pct_x10_max;   /* max |Vqd| / 32767, % x10 */
  uint16_t imag_max_c1;      /* max |I| in case 1 / case 2-3: if the spikes are only in case 2-3 they are a sampling artifact */
  uint16_t imag_max_c23;
  uint16_t imag_t_ms;        /* ms from run start (capped 65535) at the |I| maximum */
  int16_t  imag_rpm;
  uint16_t i_mag_max_s16;
  uint16_t id_abs_max_s16;
  uint16_t iq_abs_max_s16;
  int16_t  speed_at_case23_rpm;
  int16_t  speed_min_rpm;    /* negative = shaft turned backwards */
  int16_t  speed_max_rpm;
  uint16_t n_neg_ms;         /* ms with speed < -MYDIAG_NEG_RPM */
  int16_t  vd_mean_s16;      /* mean Vd, light load and speed > MYDIAG_ANG_RPM_MIN */
  int16_t  vq_mean_s16;      /* mean Vq, same condition */
  int16_t  ang_err_x10;      /* atan2(Vd,Vq) in degrees x10; ~0 = frame aligned */
  uint16_t n_ang_ms;         /* samples used */
  int16_t  rpm_mean;         /* mean speed in the window, rpm */
  int16_t  vbus_mean_x10;    /* mean bus voltage in the window, V x10 */
  int16_t  iq_std_s16;       /* ripple (standard deviation) of Iq / Id in the window, digits */
  int16_t  id_std_s16;
  int16_t  iq_mean_s16;      /* mean Iq / Id in the same window (digits) */
  int16_t  id_mean_s16;
  int16_t  pwr_mean_w_x10;   /* mean electrical power [W]x10 in the same window (compare between phase shifts) */
  int16_t  idref_mean_s16;   /* mean Id REFERENCE in the same window: if it is not 0, Id mean is not an error but a command */
} MyDiagCycle_t;

#define MYDIAG_HIST 8u
volatile MyDiagCycle_t g_myDiagHist[MYDIAG_HIST];
volatile uint32_t g_myDiagHistCount;   /* total number of finished cycles */

/* =====================================================================================
 *  PASSIVE IDENTIFICATION of R, L, lambda (Ke) and the frame angle error delta
 *
 *  Nothing is commanded, the motor is only observed. Every 64 ms of RUN the mean of
 *  Vd, Vq, Id, Iq and the speed is taken (a "window"); windows with voltage clipping,
 *  moved sampling point (case 2/3), a speed change > 8 % or no steady speed are dropped.
 *  The surface PM motor in the controller frame (isotropic, so the formula is exact
 *  for any frame offset delta):
 *
 *      Vd = R*Id - w*L(|I|)*Iq - w*lam*sin(delta)            + vdt*Id/|I|
 *      Vq = R*Iq + w*L(|I|)*Id + w*lam*cos(delta)            + vdt*Iq/|I|
 *      delta = delta0 + tau*w   (Hall interpolation + PWM delay, grows with speed)
 *      L(|I|) = L0 + L1*|I| + L2*|I|^2   (the inductance falls with current, fit[2], fit[3])
 *
 *  Unknowns (linear least squares over all windows of a phase, 2 equations per window):
 *      R [Ohm], L0 [mH], a = lam*sin(delta0), b = lam*cos(delta0) [mWb],
 *      c = lam*tau, vdt = dead-time / switch drop vector [V], L1 [mH/A], L2 [mH/A^2]
 *  Four variants are solved:
 *      fit[0] = R, L, a, b                       (simplest, R includes the dead-time drop, delta is biased by tau)
 *      fit[1] = R, L, a, b, vdt, tau             (tau only if the speed range of the windows >= 600 rpm)
 *      fit[2] = R, L0, L1, a, b, vdt, tau        (+ linear fall of L with current)
 *      fit[3] = R, L0, L1, L2, a, b, vdt, tau    (+ curvature; needs windows from ~10 A up to 100 A, otherwise noisy)
 *  Take the variant with the smallest residuals res_d/res_q whose numbers are plausible. L at 20 A and 80 A
 *  of the fitted law are in L20/L80 (equal to L for fit[0], fit[1]).
 *  All values are for the equivalent star PHASE (L = terminal value / 2, R likewise).
 *  Ke is given in the units of MOTOR_VOLTAGE_CONSTANT (V rms line-line per krpm).
 *
 *  Results (after every finished cycle, ~5 ms later; the same index as g_myDiagHist):
 *     g_myDiagIdHist[(g_myDiagHistCount-1) % 8]   or   g_myDiagIdLast
 *     g_myDiagIdAll[4]   - the same fits over ALL cycles since MyDiag_IdClearAll()
 *  delta_x100 > 0  ->  HALL_PHASE_SHIFT must be INCREASED by that angle (deg) at the mean speed
 *  of the windows (delta_ref_x100: extrapolated to MYDIAG_ID_REF_RPM, only meaningful for the fits with tau).
 *  Valid fit has n_win >= 12. If n_win == 0 the fit was not possible (too few windows or singular).
 *  IMPORTANT: the scale constants (MYDIAG_DIGITS_PER_A, MYDIAG_VQD_FS_DIV) are not verified on the
 *  real board. Compare ke_x100 with Ke measured by hand (drill + oscilloscope) once.
 * ===================================================================================== */
#if MYDIAG_ID_ENABLE

typedef struct { int16_t vd, vq, id, iq, rpm, vbus_x10; } MyDiagWin_t;

typedef struct
{
  double   ata[8][8];      /* normal equations, upper triangle used. Unknowns: 0 R, 1 L0, 2 a, 3 b, 4 c(tau), 5 vdt, 6 L1, 7 L2 */
  double   atb[8];
  double   btb;
  double   sum_om;         /* electrical speed [krad/s], summed over windows */
  double   sum_id;         /* A */
  uint32_t nw;
  int32_t  rpm_min;
  int32_t  rpm_max;
  int32_t  iq_min_x10;     /* A x10 */
  int32_t  iq_max_x10;
} MyDiagNE_t;

typedef struct
{
  int16_t n_win;           /* windows used, 0 = fit not valid */
  int16_t R_mohm_x10;      /* phase R (effective, includes wires and switches), mOhm x10 */
  int16_t L_uH_x10;        /* phase L (fit[0], fit[1]: constant; fit[2], fit[3]: L0 at zero current), uH x10 */
  int16_t L_se_uH_x10;     /* its 1-sigma, uH x10 */
  int16_t L20_uH_x10;      /* L at |I| = 20 A of the fitted law, uH x10 */
  int16_t L80_uH_x10;      /* L at |I| = 80 A */
  int16_t ke_x100;         /* V rms line-line per krpm x100 (as MOTOR_VOLTAGE_CONSTANT) */
  int16_t delta_x100;      /* frame angle error at the mean speed, deg x100; >0: increase HALL_PHASE_SHIFT */
  int16_t delta_ref_x100;  /* same at MYDIAG_ID_REF_RPM (uses tau; fits 1..3 only) */
  int16_t delta_se_x100;   /* 1-sigma of delta, deg x100 */
  int16_t tau_us;          /* speed dependent lag, us (fits 1..3) */
  int16_t vdt_mv;          /* dead-time / drop voltage, mV (fits 1..3) */
  int16_t res_d_mv;        /* rms residual of the d equation over the stored windows, mV */
  int16_t res_q_mv;        /* same for q */
} MyDiagFit_t;

typedef struct
{
  MyDiagFit_t fit[4];
  uint16_t ready;          /* 1 when finished */
  uint16_t n_win;          /* accepted windows of the phase */
  uint16_t n_rej;          /* rejected windows (clipping, ramp, case 2/3 ...) */
  uint16_t n_lost;         /* windows skipped because the previous result was still being computed */
  int16_t  rpm_min, rpm_max, rpm_mean;
  int16_t  iq_min_x10;     /* range of Iq over the windows, A x10 (the wider the better for L) */
  int16_t  iq_max_x10;
  int16_t  id_mean_x10;    /* A x10 */
  int16_t  shift_now_x10;  /* HALL_M1.PhaseShift at the moment of the report, deg x10 */
  int16_t  L_bin_uH_x10[4];/* L from the windows with |Iq| in 5-20, 20-40, 40-70, >70 A (saturation), uH x10 */
  uint16_t n_bin[4];
} MyDiagId_t;

MyDiagId_t  g_myDiagIdHist[MYDIAG_HIST];
MyDiagId_t  g_myDiagIdLast;
MyDiagFit_t g_myDiagIdAll[4];
volatile uint32_t g_myDiagIdSeq;     /* +1 every time a result is finished */

static const uint8_t s_mask[4] = { 0x0Fu, 0x3Fu, 0x7Fu, 0xFFu };   /* bit0 R, bit1 L0, bit2 a, bit3 b, bit4 c(tau), bit5 vdt, bit6 L1, bit7 L2 */
static const float   s_bin_edge[5] = { 5.0f, 20.0f, 40.0f, 70.0f, 1.0e6f };

static MyDiagNE_t  s_ph;                       /* current phase (since the last finished cycle) */
static MyDiagNE_t  s_all;                      /* all phases */
static MyDiagWin_t s_win[MYDIAG_ID_NWIN];      /* ring of the windows of the current phase */
static uint16_t    s_win_n;
static uint16_t    s_win_w;
static uint16_t    s_rej, s_lost;
static uint8_t     s_id_stage;                 /* 0 = idle, 1..11 = computing the result of the finished cycle */
static uint8_t     s_id_k;
static float       s_par[4][8];
static float       s_cov[4][8];
static bool        s_ok[4];
static uint8_t     s_msk[4];                   /* mask really used (tau dropped when the speed range is narrow) */
static MyDiagId_t  s_idout;

/* running window */
static struct
{
  int32_t  sum_vd, sum_vq, sum_id, sum_iq, sum_rpm;
  uint32_t sum_vbus;
  uint32_t n;
  int32_t  rpm_min, rpm_max;
  uint32_t vs_max;
  uint32_t dev_max;
  bool     bad;
} s_w;

static int16_t MyDiag_Sat16(float v)
{
  if (!(v == v)) { return 0; }
  if (v > 32767.0f) { return 32767; }
  if (v < -32768.0f) { return -32768; }
  return (int16_t)(v + ((v >= 0.0f) ? 0.5f : -0.5f));
}

/* window means -> physical units (V, A, krad/s electrical) and the dead-time direction */
static void MyDiag_IdPhys(const MyDiagWin_t *w, float *Vd, float *Vq, float *Id, float *Iq, float *om, float *ud, float *uq)
{
  const float vsc = (((float)w->vbus_x10 * 0.1f) / MYDIAG_VQD_FS_DIV) / 32767.0f;
  *Vd = (float)w->vd * vsc;
  *Vq = (float)w->vq * vsc;
  *Id = (float)w->id / MYDIAG_DIGITS_PER_A;
  *Iq = (float)w->iq / MYDIAG_DIGITS_PER_A;
  *om = (float)w->rpm * (6.28318531f / 60.0f) * (float)MYDIAG_POLE_PAIRS * 0.001f;
  const float im = sqrtf((*Id) * (*Id) + (*Iq) * (*Iq));
  *ud = (im > 3.0f) ? (*Id / im) : 0.0f;
  *uq = (im > 3.0f) ? (*Iq / im) : 0.0f;
}

static void MyDiag_NeAddRow(MyDiagNE_t *ne, const double r[8], double y)
{
  for (int i = 0; i < 8; i++)
  {
    ne->atb[i] += r[i] * y;
    for (int j = i; j < 8; j++)
    {
      ne->ata[i][j] += r[i] * r[j];
    }
  }
  ne->btb += y * y;
}

static void MyDiag_NeAddWindow(MyDiagNE_t *ne, const MyDiagWin_t *w)
{
  float Vd, Vq, Id, Iq, om, ud, uq;
  MyDiag_IdPhys(w, &Vd, &Vq, &Id, &Iq, &om, &ud, &uq);
  const double im = sqrt((double)Id * Id + (double)Iq * Iq);
  const double rd[8] = { Id, -(double)om * Iq, -(double)om, 0.0, -(double)om * om, ud, -(double)om * im * Iq, -(double)om * im * im * Iq };   /* unknowns: R L0 a b c vdt L1 L2 */
  const double rq[8] = { Iq,  (double)om * Id,  0.0, (double)om, 0.0, uq,  (double)om * im * Id,  (double)om * im * im * Id };
  MyDiag_NeAddRow(ne, rd, Vd);
  MyDiag_NeAddRow(ne, rq, Vq);

  const int32_t iq10 = (int32_t)(Iq * 10.0f);
  if (ne->nw == 0u)
  {
    ne->rpm_min = w->rpm; ne->rpm_max = w->rpm; ne->iq_min_x10 = iq10; ne->iq_max_x10 = iq10;
  }
  else
  {
    if (w->rpm < ne->rpm_min) { ne->rpm_min = w->rpm; }
    if (w->rpm > ne->rpm_max) { ne->rpm_max = w->rpm; }
    if (iq10 < ne->iq_min_x10) { ne->iq_min_x10 = iq10; }
    if (iq10 > ne->iq_max_x10) { ne->iq_max_x10 = iq10; }
  }
  ne->sum_om += om;
  ne->sum_id += Id;
  ne->nw++;
}

static void MyDiag_NeMerge(MyDiagNE_t *dst, const MyDiagNE_t *src)
{
  if (src->nw == 0u) { return; }
  for (int i = 0; i < 8; i++)
  {
    dst->atb[i] += src->atb[i];
    for (int j = i; j < 8; j++) { dst->ata[i][j] += src->ata[i][j]; }
  }
  dst->btb += src->btb;
  dst->sum_om += src->sum_om;
  dst->sum_id += src->sum_id;
  if (dst->nw == 0u)
  {
    dst->rpm_min = src->rpm_min; dst->rpm_max = src->rpm_max;
    dst->iq_min_x10 = src->iq_min_x10; dst->iq_max_x10 = src->iq_max_x10;
  }
  else
  {
    if (src->rpm_min < dst->rpm_min) { dst->rpm_min = src->rpm_min; }
    if (src->rpm_max > dst->rpm_max) { dst->rpm_max = src->rpm_max; }
    if (src->iq_min_x10 < dst->iq_min_x10) { dst->iq_min_x10 = src->iq_min_x10; }
    if (src->iq_max_x10 > dst->iq_max_x10) { dst->iq_max_x10 = src->iq_max_x10; }
  }
  dst->nw += src->nw;
}

/* Least squares for the unknowns selected by mask. The normal matrix is scaled to a unit diagonal and
 * inverted with Gauss-Jordan in float (FPU, microseconds). p = solution, cov = diagonal of the inverse
 * (multiply by the residual variance to get the variance of the parameter). */
static bool MyDiag_IdSolve(const MyDiagNE_t *ne, uint8_t mask, float p[8], float cov[8])
{
  uint8_t idx[8];
  uint8_t k = 0u;
  for (uint8_t i = 0u; i < 8u; i++)
  {
    p[i] = 0.0f;
    cov[i] = 0.0f;
    if (((mask >> i) & 1u) != 0u) { idx[k++] = i; }
  }
  if (ne->nw < MYDIAG_ID_MIN_WIN) { return false; }

  double d[8];
  for (uint8_t j = 0u; j < k; j++)
  {
    const double dj = ne->ata[idx[j]][idx[j]];
    if (dj < 1.0e-9) { return false; }
    d[j] = 1.0 / sqrt(dj);
  }
  float m[8][16];
  float rhs[8];
  for (uint8_t j = 0u; j < k; j++)
  {
    for (uint8_t l = 0u; l < k; l++)
    {
      const double a = (idx[j] <= idx[l]) ? ne->ata[idx[j]][idx[l]] : ne->ata[idx[l]][idx[j]];
      m[j][l] = (float)(a * d[j] * d[l]);
      m[j][k + l] = (j == l) ? 1.0f : 0.0f;
    }
    rhs[j] = (float)(ne->atb[idx[j]] * d[j]);
  }
  for (uint8_t c = 0u; c < k; c++)
  {
    uint8_t pr = c;
    for (uint8_t r = (uint8_t)(c + 1u); r < k; r++)
    {
      if (fabsf(m[r][c]) > fabsf(m[pr][c])) { pr = r; }
    }
    if (fabsf(m[pr][c]) < 1.0e-6f) { return false; }          /* singular: not enough excitation */
    if (pr != c)
    {
      for (uint8_t l = 0u; l < (uint8_t)(2u * k); l++) { const float t = m[c][l]; m[c][l] = m[pr][l]; m[pr][l] = t; }
      const float t = rhs[c]; rhs[c] = rhs[pr]; rhs[pr] = t;
    }
    const float pv = m[c][c];
    for (uint8_t l = 0u; l < (uint8_t)(2u * k); l++) { m[c][l] /= pv; }
    rhs[c] /= pv;
    for (uint8_t r = 0u; r < k; r++)
    {
      if (r == c) { continue; }
      const float f = m[r][c];
      if (f == 0.0f) { continue; }
      for (uint8_t l = 0u; l < (uint8_t)(2u * k); l++) { m[r][l] -= f * m[c][l]; }
      rhs[r] -= f * rhs[c];
    }
  }
  for (uint8_t j = 0u; j < k; j++)
  {
    p[idx[j]] = (float)((double)rhs[j] * d[j]);
    cov[idx[j]] = (float)((double)m[j][k + j] * d[j] * d[j]);
  }
  return true;
}

/* rms residual of the model over the stored windows, in volts (d and q separately) */
static float MyDiag_IdResidual(const float p[8], float *res_d, float *res_q)
{
  double sd = 0.0, sq = 0.0;
  for (uint16_t i = 0u; i < s_win_n; i++)
  {
    float Vd, Vq, Id, Iq, om, ud, uq;
    MyDiag_IdPhys(&s_win[i], &Vd, &Vq, &Id, &Iq, &om, &ud, &uq);
    const float im = sqrtf(Id * Id + Iq * Iq);
    const float L = p[1] + p[6] * im + p[7] * im * im;
    const float pd = p[0] * Id - om * L * Iq - om * p[2] - om * om * p[4] + p[5] * ud;
    const float pq = p[0] * Iq + om * L * Id + om * p[3] + p[5] * uq;
    sd += (double)((Vd - pd) * (Vd - pd));
    sq += (double)((Vq - pq) * (Vq - pq));
  }
  const float n = (float)((s_win_n > 0u) ? s_win_n : 1u);
  *res_d = sqrtf((float)sd / n);
  *res_q = sqrtf((float)sq / n);
  return ((float)(sd + sq)) / (2.0f * n);        /* residual variance per equation, V^2 */
}

static void MyDiag_IdFillFit(MyDiagFit_t *f, const MyDiagNE_t *ne, const float p[8], const float cov[8], float sigma2,
                             uint8_t mask, float res_d, float res_q)
{
  const float rad2deg = 57.2957795f;
  const float om_mean = (float)(ne->sum_om / (double)ne->nw);
  const float om_ref = (float)MYDIAG_ID_REF_RPM * (6.28318531f / 60.0f) * (float)MYDIAG_POLE_PAIRS * 0.001f;
  const float lam = hypotf(p[2], p[3]);                      /* mWb */
  f->n_win = (ne->nw > 32767u) ? 32767 : (int16_t)ne->nw;
  f->R_mohm_x10 = MyDiag_Sat16(p[0] * 1000.0f * 10.0f);
  f->L_uH_x10 = MyDiag_Sat16(p[1] * 1000.0f * 10.0f);
  f->L_se_uH_x10 = MyDiag_Sat16(sqrtf(sigma2 * cov[1]) * 1000.0f * 10.0f);
  f->L20_uH_x10 = MyDiag_Sat16((p[1] + p[6] * 20.0f + p[7] * 400.0f) * 1000.0f * 10.0f);
  f->L80_uH_x10 = MyDiag_Sat16((p[1] + p[6] * 80.0f + p[7] * 6400.0f) * 1000.0f * 10.0f);
  f->ke_x100 = MyDiag_Sat16(lam * 1.0e-3f * (6.28318531f / 60.0f) * 1000.0f * (float)MYDIAG_POLE_PAIRS * 1.22474487f * 100.0f);
  f->delta_x100 = MyDiag_Sat16(atan2f(p[2] + p[4] * om_mean, p[3]) * rad2deg * 100.0f);
  f->delta_ref_x100 = MyDiag_Sat16(atan2f(p[2] + p[4] * om_ref, p[3]) * rad2deg * 100.0f);
  f->delta_se_x100 = (lam > 0.0f) ? MyDiag_Sat16(sqrtf(sigma2 * cov[2]) / lam * rad2deg * 100.0f) : 0;
  f->tau_us = (((mask & 0x10u) != 0u) && (p[3] != 0.0f)) ? MyDiag_Sat16(p[4] / p[3] * 1000.0f) : 0;
  f->vdt_mv = MyDiag_Sat16(p[5] * 1000.0f);
  f->res_d_mv = MyDiag_Sat16(res_d * 1000.0f);
  f->res_q_mv = MyDiag_Sat16(res_q * 1000.0f);
}

static void MyDiag_IdClearPhase(void)
{
  memset(&s_ph, 0, sizeof(s_ph));
  s_win_n = 0u;
  s_win_w = 0u;
  s_rej = 0u;
  s_lost = 0u;
  s_w.n = 0u;
}

void MyDiag_IdClearAll(void)
{
  memset(&s_all, 0, sizeof(s_all));
  memset(g_myDiagIdAll, 0, sizeof(g_myDiagIdAll));
  MyDiag_IdClearPhase();
  s_id_stage = 0u;
}

static void MyDiag_IdWindowDone(void)
{
  const int32_t n = (int32_t)s_w.n;
  const int32_t rpm_mean = s_w.sum_rpm / n;
  const int32_t span = s_w.rpm_max - s_w.rpm_min;
  const bool ok = (!s_w.bad) && (rpm_mean >= MYDIAG_ID_RPM_MIN) &&
                  ((span * 100) <= (MYDIAG_ID_SPAN_PCT * rpm_mean)) &&
                  (s_w.dev_max <= (uint32_t)MYDIAG_ID_STEADY_RPM) &&
                  ((s_w.vs_max * 1000u) <= ((uint32_t)MYDIAG_ID_VS_MAX_X10 * 32767u));
  if (!ok)
  {
    s_rej++;
    return;
  }
  MyDiagWin_t w;
  w.vd = (int16_t)(s_w.sum_vd / n);
  w.vq = (int16_t)(s_w.sum_vq / n);
  w.id = (int16_t)(s_w.sum_id / n);
  w.iq = (int16_t)(s_w.sum_iq / n);
  w.rpm = (int16_t)rpm_mean;
  w.vbus_x10 = (int16_t)((s_w.sum_vbus * 10u) / s_w.n);
  s_win[s_win_w] = w;
  s_win_w = (uint16_t)((s_win_w + 1u) % MYDIAG_ID_NWIN);
  if (s_win_n < MYDIAG_ID_NWIN) { s_win_n++; }
  MyDiag_NeAddWindow(&s_ph, &w);
}

/* one call per ms, only in state RUN */
static void MyDiag_IdSample(int32_t vd, int32_t vq, int32_t id, int32_t iq, int32_t rpm, int32_t ref_rpm, uint32_t vs,
                            uint32_t vbus_v, bool bad)
{
  if (s_id_stage != 0u)
  {
    s_lost++;                    /* previous result still being computed (a few ms), skip */
    return;
  }
  if (s_w.n == 0u)
  {
    memset(&s_w, 0, sizeof(s_w));
    s_w.rpm_min = rpm;
    s_w.rpm_max = rpm;
  }
  s_w.sum_vd += vd;
  s_w.sum_vq += vq;
  s_w.sum_id += id;
  s_w.sum_iq += iq;
  s_w.sum_rpm += rpm;
  s_w.sum_vbus += vbus_v;
  if (rpm < s_w.rpm_min) { s_w.rpm_min = rpm; }
  if (rpm > s_w.rpm_max) { s_w.rpm_max = rpm; }
  if (vs > s_w.vs_max) { s_w.vs_max = vs; }
  const uint32_t dev = (uint32_t)abs(rpm - ref_rpm);
  if (dev > s_w.dev_max) { s_w.dev_max = dev; }
  if (bad) { s_w.bad = true; }
  s_w.n++;
  if (s_w.n >= MYDIAG_ID_WIN_MS)
  {
    MyDiag_IdWindowDone();
    s_w.n = 0u;
  }
}

/* variant v >= 1 solves for tau only when the windows cover a wide speed range (otherwise tau is not separable from delta0) */
static uint8_t MyDiag_IdMask(uint8_t v, const MyDiagNE_t *ne)
{
  uint8_t m = s_mask[v];
  if ((v >= 1u) && ((ne->rpm_max - ne->rpm_min) < MYDIAG_ID_TAU_SPAN_RPM)) { m &= (uint8_t)~0x10u; }
  return m;
}

/* One stage per call (each < 0.2 ms) so that the 1 kHz task never overruns. Called from Step() in all
 * states except RUN and from the early-return path of CycleEnd(), i.e. many times in the Stop state. */
void MyDiag_IdWork(void)
{
  switch (s_id_stage)
  {
    case 0:
      return;

    case 1:   /* statistics of the phase that has just finished */
    {
      memset(&s_idout, 0, sizeof(s_idout));
      memset(s_ok, 0, sizeof(s_ok));
      s_idout.n_win = (s_ph.nw > 65535u) ? 65535u : (uint16_t)s_ph.nw;
      s_idout.n_rej = s_rej;
      s_idout.n_lost = s_lost;
      if (s_ph.nw > 0u)
      {
        s_idout.rpm_min = (int16_t)s_ph.rpm_min;
        s_idout.rpm_max = (int16_t)s_ph.rpm_max;
        s_idout.rpm_mean = MyDiag_Sat16((float)(s_ph.sum_om * 1000.0 / (6.28318531 / 60.0 * (double)MYDIAG_POLE_PAIRS) / (double)s_ph.nw));
        s_idout.iq_min_x10 = MyDiag_Sat16((float)s_ph.iq_min_x10);
        s_idout.iq_max_x10 = MyDiag_Sat16((float)s_ph.iq_max_x10);
        s_idout.id_mean_x10 = MyDiag_Sat16((float)(s_ph.sum_id / (double)s_ph.nw * 10.0));
      }
#if MYDIAG_USE_HALL_M1
      {
        int32_t sh = ((int32_t)HALL_M1.PhaseShift * 3600) / 65536;
        if (sh < 0) { sh += 3600; }
        s_idout.shift_now_x10 = (int16_t)sh;
      }
#endif
      s_id_stage = (s_ph.nw >= MYDIAG_ID_MIN_WIN) ? 2u : 6u;
      break;
    }

    case 2:
    case 3:
    case 4:
    case 5:   /* variant 0, 1, 2, 3: solve, residual (sigma), final numbers */
    {
      const uint8_t v = (uint8_t)(s_id_stage - 2u);
      s_msk[v] = MyDiag_IdMask(v, &s_ph);
      s_ok[v] = MyDiag_IdSolve(&s_ph, s_msk[v], s_par[v], s_cov[v]);
      if (s_ok[v])
      {
        float rd, rq;
        const float sigma2 = MyDiag_IdResidual(s_par[v], &rd, &rq);
        MyDiag_IdFillFit(&s_idout.fit[v], &s_ph, s_par[v], s_cov[v], sigma2, s_msk[v], rd, rq);
      }
      s_id_stage++;
      break;
    }

    case 6:   /* saturation bins of the best variant, add the phase to the all-cycles sums */
    {
      const int8_t bv = s_ok[3] ? 3 : (s_ok[2] ? 2 : (s_ok[1] ? 1 : (s_ok[0] ? 0 : -1)));
      if (bv >= 0)
      {
        const float *p = s_par[bv];
        float sum[4] = { 0.0f, 0.0f, 0.0f, 0.0f };
        uint16_t cnt[4] = { 0u, 0u, 0u, 0u };
        for (uint16_t i = 0u; i < s_win_n; i++)
        {
          float Vd, Vq, Id, Iq, om, ud, uq;
          MyDiag_IdPhys(&s_win[i], &Vd, &Vq, &Id, &Iq, &om, &ud, &uq);
          const float aiq = fabsf(Iq);
          if ((aiq < s_bin_edge[0]) || (om < 0.2f)) { continue; }
          const float lmh = (p[0] * Id - om * p[2] - om * om * p[4] + p[5] * ud - Vd) / (om * Iq);   /* apparent L of this window */
          for (uint8_t b = 0u; b < 4u; b++)
          {
            if ((aiq >= s_bin_edge[b]) && (aiq < s_bin_edge[b + 1u])) { sum[b] += lmh; cnt[b]++; break; }
          }
        }
        for (uint8_t b = 0u; b < 4u; b++)
        {
          s_idout.n_bin[b] = cnt[b];
          s_idout.L_bin_uH_x10[b] = (cnt[b] > 0u) ? MyDiag_Sat16(sum[b] / (float)cnt[b] * 1.0e4f) : 0;
        }
      }
      MyDiag_NeMerge(&s_all, &s_ph);
      s_id_stage = 7u;
      break;
    }

    case 7:
    case 8:
    case 9:
    case 10:  /* the same four fits over all phases (no sigma / residual: those fields stay 0) */
    {
      const uint8_t v = (uint8_t)(s_id_stage - 7u);
      float p[8], cv[8];
      memset(&g_myDiagIdAll[v], 0, sizeof(MyDiagFit_t));
      const uint8_t mk = MyDiag_IdMask(v, &s_all);
      if (MyDiag_IdSolve(&s_all, mk, p, cv)) { MyDiag_IdFillFit(&g_myDiagIdAll[v], &s_all, p, cv, 0.0f, mk, 0.0f, 0.0f); }
      s_id_stage++;
      break;
    }

    case 11:
    default:  /* publish */
    {
      s_idout.ready = 1u;
      g_myDiagIdHist[s_id_k] = s_idout;
      g_myDiagIdLast = s_idout;
      g_myDiagIdSeq++;
      MyDiag_IdClearPhase();
      s_id_stage = 0u;
      break;
    }
  }
}

/* called by CycleEnd() when a phase has been recorded in g_myDiagHist[k] */
static void MyDiag_IdBegin(uint8_t k)
{
  if (s_id_stage != 0u) { return; }     /* previous result not finished: this phase is not evaluated */
  s_id_k = k;
  g_myDiagIdHist[k].ready = 0u;
  s_id_stage = 1u;
}

#define MyDiag_IdIdle()  do { s_w.n = 0u; MyDiag_IdWork(); } while (0)

#else  /* !MYDIAG_ID_ENABLE */
#define MyDiag_IdIdle()                          ((void)0)
#define MyDiag_IdBegin(k)                        ((void)(k))
#define MyDiag_IdSample(a, b, c, d, e, f, g, h, i) ((void)0)
#define MyDiag_IdClearAll()                      ((void)0)
#define MyDiag_IdWork()                          ((void)0)
#endif


/* =====================================================================================
 *  HALL SECTOR DIAGNOSTICS
 *
 *  With HALL_MTPA = true the library (hall_speed_pos_fdbk.c) does this:
 *    - at EVERY Hall edge hElAngle is overwritten with the nominal edge angle PhaseShift + k*60 deg;
 *    - between the edges the angle is integrated with the average speed = 60 deg / (mean of the last
 *      HALL_AVERAGING_FIFO_DEPTH sector periods) and is clamped to +60 deg until the next edge.
 *  So the angle error is (a) the placement error of the edge that starts the sector (the sensors / magnet are
 *  not exactly 60 deg apart) and (b) a ramp inside the sector, which appears when the averaging depth is not a
 *  multiple of 6 (the mean of 4 or 16 sector periods is not the mean of a whole revolution).
 *
 *  The module reads HALL_M1.SensorPeriod[] (the library's own FIFO) from the 1 kHz task, takes only NEW entries,
 *  labels them with the Hall state, and for every complete revolution at steady speed computes the width of each
 *  sector in electrical degrees (nominal 60.00). From the mean widths it
 *    - gives the edge errors edge_x100[6] (the table that would remove (a)),
 *    - PREDICTS the rms and peak-to-peak angle error that the library makes with the present averaging depth and with 6, 12, 18.
 *  Nothing is changed in the library. Result per phase: g_myDiagHallHist[same index as g_myDiagHist], g_myDiagHallLast.
 *  Sectors are in the state order 5,1,3,2,6,4 (forward rotation); edge_x100[i] is the error of the edge that
 *  STARTS sector i = true angle - nominal angle, deg x100, mean removed.
 *  Compare pred_rms_x100[0] with the measured ripple: id_std / iq_mean = sin(error). If the prediction is much smaller
 *  than the measured ripple, the Hall edges are not the main cause.
 * ===================================================================================== */
#if MYDIAG_HALL_ENABLE
#if !MYDIAG_ID_ENABLE
#error MYDIAG_HALL_ENABLE needs MYDIAG_ID_ENABLE (MyDiag_Sat16)
#endif

typedef struct
{
  uint16_t ready;
  uint16_t n_rev;               /* revolutions used */
  uint16_t n_skip;              /* revolutions dropped (speed changed) */
  uint16_t n_lost;              /* sync losses (missed edges / direction change) */
  uint8_t  N;                   /* HALL_M1.SpeedBufferSize = HALL_AVERAGING_FIFO_DEPTH */
  uint8_t  pad;
  int16_t  rpm_mean;
  int16_t  w_mean_x100[6];      /* sector widths, deg x100, nominal 6000 */
  int16_t  w_std_x100[6];       /* revolution to revolution spread of the width */
  int16_t  edge_x100[6];        /* edge errors (true - nominal), deg x100, mean removed */
  int16_t  edge_rms_x100;       /* rms of the edge errors */
  int16_t  pred_rms_x100[4];    /* predicted rms angle error for the depth N, 6, 12, 18 (deg x100) */
  int16_t  pred_pp_x100[4];     /* predicted peak to peak */
} MyDiagHall_t;

MyDiagHall_t g_myDiagHallHist[MYDIAG_HIST];
MyDiagHall_t g_myDiagHallLast;

static const uint8_t s_hl_pos[8] = { 255u, 1u, 3u, 2u, 5u, 0u, 4u, 255u };   /* state -> position in the sequence 5,1,3,2,6,4 */

static uint8_t  s_hl_last_idx;
static bool     s_hl_have;
static uint8_t  s_hl_last_dpos;
static uint32_t s_hl_cnt;
static uint32_t s_hl_d[6];
static uint32_t s_hl_Sprev;
static float    s_hl_sum[6];
static float    s_hl_sum2[6];
static uint32_t s_hl_nrev, s_hl_nskip, s_hl_nlost;
static float    s_hl_sum_rpm;
static int32_t  s_hl_rpm;

static void MyDiag_HallClear(void)
{
  s_hl_have = false;
  s_hl_cnt = 0u;
  s_hl_Sprev = 0u;
  memset(s_hl_sum, 0, sizeof(s_hl_sum));
  memset(s_hl_sum2, 0, sizeof(s_hl_sum2));
  s_hl_nrev = 0u;
  s_hl_nskip = 0u;
  s_hl_nlost = 0u;
  s_hl_sum_rpm = 0.0f;
}

static void MyDiag_HallRev(void)
{
  uint32_t S = 0u;
  for (uint8_t p = 0u; p < 6u; p++) { S += s_hl_d[p]; }
  if (S == 0u) { return; }
  const uint32_t Sp = s_hl_Sprev;
  s_hl_Sprev = S;
  if (Sp != 0u)
  {
    const uint32_t diff = (S > Sp) ? (S - Sp) : (Sp - S);
    if ((diff * 100u) > ((uint32_t)MYDIAG_HALL_REV_TOL_PCT * S)) { s_hl_nskip++; return; }
  }
  else
  {
    return;                       /* first revolution after a (re)start: no reference yet */
  }
  for (uint8_t p = 0u; p < 6u; p++)
  {
    const float e = 36000.0f * (float)s_hl_d[p] / (float)S - 6000.0f;      /* deg x100 */
    s_hl_sum[p] += e;
    s_hl_sum2[p] += e * e;
  }
  s_hl_nrev++;
  s_hl_sum_rpm += (float)s_hl_rpm;
}

/* called every ms in RUN */
static void MyDiag_HallStep(int32_t rpm, int32_t ref_rpm)
{
  const uint8_t N = HALL_M1.SpeedBufferSize;
  if ((N < 2u) || (N > 18u)) { return; }
  const uint8_t idx1 = *(volatile uint8_t *)&HALL_M1.SpeedFIFOIdx;
  const uint8_t st1 = *(volatile uint8_t *)&HALL_M1.HallState;
  const bool steady = (rpm >= MYDIAG_HALL_RPM_MIN) && (abs(rpm - ref_rpm) < MYDIAG_HALL_STEADY_RPM) &&
                      (*(volatile int8_t *)&HALL_M1.Direction > 0);
  if ((!steady) || (idx1 >= N))
  {
    s_hl_last_idx = idx1;
    s_hl_have = false;
    s_hl_cnt = 0u;
    s_hl_Sprev = 0u;
    return;
  }
  if (!s_hl_have)
  {
    s_hl_last_idx = idx1;
    s_hl_have = true;
    s_hl_cnt = 0u;
    s_hl_Sprev = 0u;
    return;
  }
  const uint8_t m = (uint8_t)((idx1 + N - s_hl_last_idx) % N);
  if (m == 0u) { return; }
  const uint8_t pos = (st1 < 8u) ? s_hl_pos[st1] : 255u;
  if (pos == 255u)
  {
    s_hl_last_idx = idx1;
    s_hl_cnt = 0u;
    s_hl_Sprev = 0u;
    return;
  }
  int32_t buf[18];
  for (uint8_t a = 0u; a < m; a++)
  {
    const uint8_t i = (uint8_t)((idx1 + N - 1u - a) % N);
    int32_t v = ((volatile int32_t *)HALL_M1.SensorPeriod)[i];
    buf[a] = (v < 0) ? -v : v;
  }
  if ((*(volatile uint8_t *)&HALL_M1.SpeedFIFOIdx != idx1) || (*(volatile uint8_t *)&HALL_M1.HallState != st1))
  {
    return;                       /* an edge came while we were reading: try again in 1 ms (m will be larger) */
  }
  s_hl_last_idx = idx1;
  s_hl_rpm = rpm;
  for (int32_t a = (int32_t)m - 1; a >= 0; a--)    /* oldest new entry first */
  {
    /* entry of age a was stored when the state at position (pos-a) was entered; its value is the duration of the state before it */
    const uint8_t dpos = (uint8_t)((pos + 5u + 6u - ((uint32_t)a % 6u)) % 6u);
    if ((s_hl_cnt > 0u) && (dpos != (uint8_t)((s_hl_last_dpos + 1u) % 6u)))
    {
      s_hl_nlost++;
      s_hl_cnt = 0u;
      s_hl_Sprev = 0u;
    }
    s_hl_d[dpos] = (uint32_t)buf[a];
    s_hl_last_dpos = dpos;
    s_hl_cnt++;
    if ((s_hl_cnt >= 6u) && ((s_hl_cnt % 6u) == 0u)) { MyDiag_HallRev(); }
  }
}

/* predicted angle error of the library for the averaging depth N, from the mean sector widths (deg) and edge errors (deg) */
static void MyDiag_HallPredict(const float width[6], const float edge[6], uint8_t N, float *rms, float *pp)
{
  float sw = 0.0f, sm = 0.0f, mn = 1.0e9f, mx = -1.0e9f;
  float err[72];
  float wt[6];
  for (uint8_t p = 0u; p < 6u; p++)
  {
    float wb = 0.0f;
    for (uint8_t k = 1u; k <= N; k++) { wb += width[(p + 36u - k) % 6u]; }
    wb /= (float)N;
    wt[p] = width[p];
    for (uint8_t j = 0u; j < 12u; j++)
    {
      const float sa = ((float)j + 0.5f) / 12.0f * width[p];       /* true angle travelled in the sector */
      float adv = sa * 60.0f / wb;                                 /* library: average speed integration ... */
      if (adv > 60.0f) { adv = 60.0f; }                            /* ... clamped at the next edge angle */
      const float e = adv - sa - edge[p];
      err[p * 12u + j] = e;
      sm += e * width[p];
      sw += width[p];
      if (e < mn) { mn = e; }
      if (e > mx) { mx = e; }
    }
  }
  const float mean = sm / sw;
  float sq = 0.0f;
  for (uint8_t p = 0u; p < 6u; p++)
  {
    for (uint8_t j = 0u; j < 12u; j++)
    {
      const float d = err[p * 12u + j] - mean;
      sq += d * d * wt[p];
    }
  }
  *rms = sqrtf(sq / sw);
  *pp = mx - mn;
}

static void MyDiag_HallFinish(uint8_t k)
{
  MyDiagHall_t o;
  memset(&o, 0, sizeof(o));
  o.N = HALL_M1.SpeedBufferSize;
  o.n_rev = (s_hl_nrev > 65535u) ? 65535u : (uint16_t)s_hl_nrev;
  o.n_skip = (s_hl_nskip > 65535u) ? 65535u : (uint16_t)s_hl_nskip;
  o.n_lost = (s_hl_nlost > 65535u) ? 65535u : (uint16_t)s_hl_nlost;
  if (s_hl_nrev >= 10u)
  {
    const float n = (float)s_hl_nrev;
    float width[6], edge[6], e[6];
    float cum = 0.0f, cmean = 0.0f;
    o.rpm_mean = MyDiag_Sat16(s_hl_sum_rpm / n);
    for (uint8_t p = 0u; p < 6u; p++)
    {
      const float m = s_hl_sum[p] / n;
      float v = s_hl_sum2[p] / n - m * m;
      if (v < 0.0f) { v = 0.0f; }
      o.w_mean_x100[p] = MyDiag_Sat16(6000.0f + m);
      o.w_std_x100[p] = MyDiag_Sat16(sqrtf(v));
      width[p] = 60.0f + m / 100.0f;
      e[p] = cum;                         /* E_p = sum of the deviations of the sectors before p, deg */
      cmean += cum;
      cum += m / 100.0f;
    }
    cmean /= 6.0f;
    float sq = 0.0f;
    for (uint8_t p = 0u; p < 6u; p++)
    {
      edge[p] = e[p] - cmean;
      o.edge_x100[p] = MyDiag_Sat16(edge[p] * 100.0f);
      sq += edge[p] * edge[p];
    }
    o.edge_rms_x100 = MyDiag_Sat16(sqrtf(sq / 6.0f) * 100.0f);
    const uint8_t ns[4] = { (o.N >= 2u) ? o.N : 6u, 6u, 12u, 18u };
    for (uint8_t i = 0u; i < 4u; i++)
    {
      float r, pp;
      MyDiag_HallPredict(width, edge, ns[i], &r, &pp);
      o.pred_rms_x100[i] = MyDiag_Sat16(r * 100.0f);
      o.pred_pp_x100[i] = MyDiag_Sat16(pp * 100.0f);
    }
  }
  o.ready = 1u;
  g_myDiagHallHist[k] = o;
  g_myDiagHallLast = o;
  MyDiag_HallClear();
}

#define MyDiag_HallIdle()   do { s_hl_have = false; s_hl_cnt = 0u; s_hl_Sprev = 0u; } while (0)

#else  /* !MYDIAG_HALL_ENABLE */
#define MyDiag_HallStep(a, b)   ((void)0)
#define MyDiag_HallFinish(k)    ((void)(k))
#define MyDiag_HallIdle()       ((void)0)
#define MyDiag_HallClear()      ((void)0)
#endif

static PWMC_R3_2_Handle_t *MyDiag_Pwm(void)
{
  return (PWMC_R3_2_Handle_t *)pwmcHandle[MYDIAG_MOTOR_IDX];  /* cstat !MISRAC2012-Rule-11.3 */
}

void MyDiag_Reset(void)
{
  PWMC_R3_2_Handle_t *h = MyDiag_Pwm();

  g_myDiag.n_total = 0u;
  g_myDiag.n_case1 = 0u;
  g_myDiag.n_case23 = 0u;
  g_myDiag.n_est = 0u;
  g_myDiag.pct_case23_x10 = 0u;
  g_myDiag.pct_est_x10 = 0u;
  g_myDiag.headroom_min = 0x7FFFFFFF;
  g_myDiag.headroom_ns_min = 0x7FFFFFFF;
  g_myDiag.lowDuty_max = 0u;
  g_myDiag.duty_pct_x10_max = 0u;
  g_myDiag.vs_s16_max = 0u;
  g_myDiag.vs_pct_x10_max = 0u;
  g_myDiag.iq_abs_max_s16 = 0u;
  g_myDiag.i_mag_max_s16 = 0u;
  g_myDiag.imag_max_c1 = 0u;
  g_myDiag.imag_max_c23 = 0u;
  g_myDiag.imag_t_ms = 0u;
  g_myDiag.imag_rpm = 0;
  g_myDiag.id_abs_max_s16 = 0u;
  g_myDiag.iq_at_case23 = 0;
  g_myDiag.vs_at_case23 = 0u;
  g_myDiag.speed_at_case23_rpm = 0;
  g_myDiag.speed_min_rpm = 32767;
  g_myDiag.speed_max_rpm = -32768;
  g_myDiag.n_neg = 0u;
  g_myDiag.sum_vd = 0;
  g_myDiag.sum_vq = 0;
  g_myDiag.n_ang = 0u;
  g_myDiag.sum_pw = 0.0f;
  g_myDiag.sum_iq = 0;
  g_myDiag.sum_rpm = 0.0f;
  g_myDiag.sum_vbus = 0.0f;
  g_myDiag.sum_iq2 = 0.0f;
  g_myDiag.sum_id2 = 0.0f;
  g_myDiag.sum_id = 0;
  g_myDiag.sum_idref = 0;

  if ((h != NULL) && (h->pParams_str != NULL))
  {
    g_myDiag.Tafter = h->pParams_str->Tafter;
    g_myDiag.Half = h->Half_PWMPeriod;
    g_myDiag.Tbefore = h->pParams_str->Tbefore;
    g_myDiag.Tcase2 = h->pParams_str->Tcase2;
    g_myDiag.Tcase3 = h->pParams_str->Tcase3;
    g_myDiag.Tsampling = h->pParams_str->Tsampling;
  }
}

void MyDiag_Init(void)
{
  MyDiag_Reset();
  MyDiag_IdClearAll();
  MyDiag_HallClear();
}

/* Call as often as possible. Takes < 1 us, no blocking. */
void MyDiag_Step(void)
{
  PWMC_R3_2_Handle_t *h = MyDiag_Pwm();
  if ((h == NULL) || (h->pParams_str == NULL))
  {
    return;
  }

  /* count only while the motor really runs (otherwise Stop/Wait pollute the statistic) */
  if (MCI_GetSTMState(pMCI[MYDIAG_MOTOR_IDX]) != RUN)
  {
    MyDiag_IdIdle();   /* drops the unfinished window and advances the result computation, < 0.2 ms */
    MyDiag_HallIdle();
    return;
  }

  const uint32_t half   = h->Half_PWMPeriod;
  const uint32_t tafter = h->pParams_str->Tafter;
  const uint32_t low    = h->_Super.lowDuty;
  const uint32_t ccr4   = LL_TIM_OC_GetCompareCH4(h->pParams_str->TIMx);

  /* headroom of case 1 condition: Half - lowDuty > Tafter */
  const int32_t headroom = (int32_t)half - (int32_t)low - (int32_t)tafter;
  const int32_t headroom_ns = (headroom * 1000) / (int32_t)MYDIAG_TIM_CLK_MHZ;

  /* real sampling case: case 1 only if CCR4 is exactly in the middle */
  const bool case1 = (ccr4 == (half - 1u));
  const bool est   = h->_Super.useEstCurrent;

  qd_t Vqd = MCI_GetVqd(pMCI[MYDIAG_MOTOR_IDX]);
  qd_t Iqd = MCI_GetIqd(pMCI[MYDIAG_MOTOR_IDX]);
  const int32_t vq = (int32_t)Vqd.q;
  const int32_t vd = (int32_t)Vqd.d;
  /* integer sqrt via float, ok for a diagnostic */
  const uint32_t vs = (uint32_t)__builtin_sqrtf((float)(vq * vq) + (float)(vd * vd));

  g_myDiag.n_total++;
  {
    const int32_t rpm = (int32_t)MCI_GetAvrgMecSpeedUnit(pMCI[MYDIAG_MOTOR_IDX]) * 6;
    if (rpm < g_myDiag.speed_min_rpm) { g_myDiag.speed_min_rpm = rpm; }
    if (rpm > g_myDiag.speed_max_rpm) { g_myDiag.speed_max_rpm = rpm; }
    if (rpm < -MYDIAG_NEG_RPM) { g_myDiag.n_neg++; }
    const int32_t ref_rpm = (int32_t)MCI_GetMecSpeedRefUnit(pMCI[MYDIAG_MOTOR_IDX]) * 6;
    if ((rpm > MYDIAG_ANG_RPM_MIN) && (abs(rpm - ref_rpm) < MYDIAG_ANG_STEADY_RPM) && (abs((int32_t)Iqd.q) < MYDIAG_ANG_IQ_MAX) &&
        (abs((int32_t)Iqd.d) < MYDIAG_ANG_IQ_MAX))
    {
      g_myDiag.sum_vd += vd;
      g_myDiag.sum_vq += vq;
      g_myDiag.n_ang++;
      g_myDiag.sum_iq += (int32_t)Iqd.q;
      g_myDiag.sum_rpm += (float)rpm;
      g_myDiag.sum_vbus += (float)VBS_GetAvBusVoltage_V(&BusVoltageSensor_M1._Super);
      g_myDiag.sum_iq2 += (float)((int32_t)Iqd.q * (int32_t)Iqd.q);
      g_myDiag.sum_id2 += (float)((int32_t)Iqd.d * (int32_t)Iqd.d);
      g_myDiag.sum_id += (int32_t)Iqd.d;
      g_myDiag.sum_pw += PQD_GetAvrgElMotorPowerW(pMPM[MYDIAG_MOTOR_IDX]);
#if MYDIAG_USE_IDREF
      g_myDiag.sum_idref += (int32_t)MCI_GetIqdref(pMCI[MYDIAG_MOTOR_IDX]).d;
#endif
    }
    MyDiag_IdSample(vd, vq, (int32_t)Iqd.d, (int32_t)Iqd.q, rpm, ref_rpm, vs,
                    (uint32_t)VBS_GetAvBusVoltage_V(&BusVoltageSensor_M1._Super), (!case1) || est);
    MyDiag_HallStep(rpm, ref_rpm);
  }
  if (case1)
  {
    g_myDiag.n_case1++;
  }
  else
  {
    g_myDiag.n_case23++;
    g_myDiag.iq_at_case23 = (int32_t)Iqd.q;
    g_myDiag.vs_at_case23 = vs;
    g_myDiag.speed_at_case23_rpm = ((int32_t)MCI_GetAvrgMecSpeedUnit(pMCI[MYDIAG_MOTOR_IDX]) * 6);
  }
  if (est)
  {
    g_myDiag.n_est++;
  }

  if (headroom < g_myDiag.headroom_min)
  {
    g_myDiag.headroom_min = headroom;
    g_myDiag.headroom_ns_min = headroom_ns;
  }
  if (low > g_myDiag.lowDuty_max)
  {
    g_myDiag.lowDuty_max = low;
    g_myDiag.duty_pct_x10_max = (half != 0u) ? ((low * 1000u) / half) : 0u;
  }
  if (vs > g_myDiag.vs_s16_max)
  {
    g_myDiag.vs_s16_max = vs;
    g_myDiag.vs_pct_x10_max = (vs * 1000u) / 32767u;
  }
  {
    const uint32_t ida = (uint32_t)abs((int32_t)Iqd.d);
    const uint32_t imag = (uint32_t)__builtin_sqrtf((float)((int32_t)Iqd.q * (int32_t)Iqd.q) +
                                                    (float)((int32_t)Iqd.d * (int32_t)Iqd.d));
    if (ida > g_myDiag.id_abs_max_s16) { g_myDiag.id_abs_max_s16 = ida; }
    if (imag > g_myDiag.i_mag_max_s16)
    {
      g_myDiag.i_mag_max_s16 = imag;
      g_myDiag.imag_t_ms = g_myDiag.n_total;
      g_myDiag.imag_rpm = (int32_t)MCI_GetAvrgMecSpeedUnit(pMCI[MYDIAG_MOTOR_IDX]) * 6;
    }
    if (case1)
    {
      if (imag > g_myDiag.imag_max_c1) { g_myDiag.imag_max_c1 = imag; }
    }
    else
    {
      if (imag > g_myDiag.imag_max_c23) { g_myDiag.imag_max_c23 = imag; }
    }
    const uint32_t iqa = (uint32_t)abs((int32_t)Iqd.q);
    if (iqa > g_myDiag.iq_abs_max_s16)
    {
      g_myDiag.iq_abs_max_s16 = iqa;
    }
  }

  g_myDiag.pct_case23_x10 = (g_myDiag.n_case23 * 1000u) / g_myDiag.n_total;
  g_myDiag.pct_est_x10 = (g_myDiag.n_est * 1000u) / g_myDiag.n_total;
}

/* Call ONCE at the end of a pump cycle (e.g. in state Stop), then it starts counting again.
 * The results of the finished cycle are in g_myDiagHist[(g_myDiagHistCount-1) % 8]. */
void MyDiag_CycleEnd(void)
{
  /* Stop state is executed ~1000 times (relief pulse): record only real runs, once */
  if (g_myDiag.n_total < MYDIAG_MIN_SAMPLES)
  {
    MyDiag_IdWork();   /* the Stop state calls this ~1000 times: the identification result is computed here */
    return;
  }
  volatile MyDiagCycle_t *c = &g_myDiagHist[g_myDiagHistCount % MYDIAG_HIST];
  int32_t hr = g_myDiag.headroom_ns_min;
  if (hr > 32767) { hr = 32767; }
  if (hr < -32768) { hr = -32768; }
  c->n_total = g_myDiag.n_total;
  c->pct_case23_x10 = (uint16_t)g_myDiag.pct_case23_x10;
  c->pct_est_x10 = (uint16_t)g_myDiag.pct_est_x10;
  c->headroom_ns_min = (int16_t)hr;
  c->duty_pct_x10_max = (uint16_t)g_myDiag.duty_pct_x10_max;
  c->vs_pct_x10_max = (uint16_t)g_myDiag.vs_pct_x10_max;
  c->iq_abs_max_s16 = (uint16_t)g_myDiag.iq_abs_max_s16;
  c->imag_max_c1 = (uint16_t)((g_myDiag.imag_max_c1 > 65535u) ? 65535u : g_myDiag.imag_max_c1);
  c->imag_max_c23 = (uint16_t)((g_myDiag.imag_max_c23 > 65535u) ? 65535u : g_myDiag.imag_max_c23);
  c->imag_t_ms = (uint16_t)((g_myDiag.imag_t_ms > 65535u) ? 65535u : g_myDiag.imag_t_ms);
  c->imag_rpm = (int16_t)g_myDiag.imag_rpm;
  c->i_mag_max_s16 = (uint16_t)((g_myDiag.i_mag_max_s16 > 65535u) ? 65535u : g_myDiag.i_mag_max_s16);
  c->id_abs_max_s16 = (uint16_t)((g_myDiag.id_abs_max_s16 > 65535u) ? 65535u : g_myDiag.id_abs_max_s16);
  c->speed_at_case23_rpm = (int16_t)g_myDiag.speed_at_case23_rpm;
  c->speed_min_rpm = (int16_t)g_myDiag.speed_min_rpm;
  c->speed_max_rpm = (int16_t)g_myDiag.speed_max_rpm;
  c->n_neg_ms = (g_myDiag.n_neg > 65535u) ? 65535u : (uint16_t)g_myDiag.n_neg;
  if (g_myDiag.n_ang >= 50u)
  {
    const float vdm = (float)g_myDiag.sum_vd / (float)g_myDiag.n_ang;
    const float vqm = (float)g_myDiag.sum_vq / (float)g_myDiag.n_ang;
    c->vd_mean_s16 = (int16_t)vdm;
    c->vq_mean_s16 = (int16_t)vqm;
    c->ang_err_x10 = (int16_t)(atan2f(vdm, vqm) * (1800.0f / 3.14159265f));
    c->pwr_mean_w_x10 = (int16_t)(g_myDiag.sum_pw / (float)g_myDiag.n_ang * 10.0f);
    c->iq_mean_s16 = (int16_t)(g_myDiag.sum_iq / (int32_t)g_myDiag.n_ang);
    c->rpm_mean = (int16_t)(g_myDiag.sum_rpm / (float)g_myDiag.n_ang);
    c->vbus_mean_x10 = (int16_t)(g_myDiag.sum_vbus / (float)g_myDiag.n_ang * 10.0f);
    c->id_mean_s16 = (int16_t)(g_myDiag.sum_id / (int32_t)g_myDiag.n_ang);
    c->idref_mean_s16 = (int16_t)(g_myDiag.sum_idref / (int32_t)g_myDiag.n_ang);
    {
      const float n = (float)g_myDiag.n_ang;
      const float mq = (float)g_myDiag.sum_iq / n;
      const float md = (float)g_myDiag.sum_id / n;
      float vq2 = g_myDiag.sum_iq2 / n - mq * mq;
      float vd2 = g_myDiag.sum_id2 / n - md * md;
      if (vq2 < 0.0f) { vq2 = 0.0f; }
      if (vd2 < 0.0f) { vd2 = 0.0f; }
      c->iq_std_s16 = (int16_t)sqrtf(vq2);
      c->id_std_s16 = (int16_t)sqrtf(vd2);
    }
  }
  else
  {
    c->vd_mean_s16 = 0;
    c->vq_mean_s16 = 0;
    c->ang_err_x10 = 0;
    c->pwr_mean_w_x10 = 0;
    c->iq_mean_s16 = 0;
    c->rpm_mean = 0;
    c->vbus_mean_x10 = 0;
    c->id_mean_s16 = 0;
    c->iq_std_s16 = 0;
    c->id_std_s16 = 0;
    c->idref_mean_s16 = 0;
  }
  c->n_ang_ms = (g_myDiag.n_ang > 65535u) ? 65535u : (uint16_t)g_myDiag.n_ang;
  MyDiag_HallFinish((uint8_t)(g_myDiagHistCount % MYDIAG_HIST));
  MyDiag_IdBegin((uint8_t)(g_myDiagHistCount % MYDIAG_HIST));   /* identification result goes to g_myDiagIdHist[same index] */
  g_myDiagHistCount++;
  MyDiag_Reset();
}

/* Live estimate of the angle error in degrees x10 (0 if not enough samples yet). */
int32_t MyDiag_AngleErrX10(void)
{
  if (g_myDiag.n_ang < 50u) { return 0; }
  const float vdm = (float)g_myDiag.sum_vd / (float)g_myDiag.n_ang;
  const float vqm = (float)g_myDiag.sum_vq / (float)g_myDiag.n_ang;
  return (int32_t)(atan2f(vdm, vqm) * (1800.0f / 3.14159265f));
}

/* Handy single numbers of the last finished cycle for a spare log channel. */
float MyDiag_LastCase23Percent(void)
{
  if (g_myDiagHistCount == 0u) { return 0.0f; }
  return (float)g_myDiagHist[(g_myDiagHistCount - 1u) % MYDIAG_HIST].pct_case23_x10 / 10.0f;
}

float MyDiag_LastHeadroomUs(void)
{
  if (g_myDiagHistCount == 0u) { return 0.0f; }
  return (float)g_myDiagHist[(g_myDiagHistCount - 1u) % MYDIAG_HIST].headroom_ns_min / 1000.0f;
}

/* One text line for a log. Returns the number of characters. */
int MyDiag_Report(char *buf, unsigned size)
{
  return snprintf(buf, size,
                  "SMP n=%lu c23=%lu.%lu%% est=%lu.%lu%% hr_min=%ld(%ldns) duty_max=%lu.%lu%% vs_max=%lu.%lu%% "
                  "Tafter=%lu Half=%lu",
                  (unsigned long)g_myDiag.n_total,
                  (unsigned long)(g_myDiag.pct_case23_x10 / 10u), (unsigned long)(g_myDiag.pct_case23_x10 % 10u),
                  (unsigned long)(g_myDiag.pct_est_x10 / 10u), (unsigned long)(g_myDiag.pct_est_x10 % 10u),
                  (long)g_myDiag.headroom_min, (long)g_myDiag.headroom_ns_min,
                  (unsigned long)(g_myDiag.duty_pct_x10_max / 10u), (unsigned long)(g_myDiag.duty_pct_x10_max % 10u),
                  (unsigned long)(g_myDiag.vs_pct_x10_max / 10u), (unsigned long)(g_myDiag.vs_pct_x10_max % 10u),
                  (unsigned long)g_myDiag.Tafter, (unsigned long)g_myDiag.Half);
}
