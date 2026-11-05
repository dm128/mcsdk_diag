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
 * TEST PLAN FOR TRISE
 *   run the same cycle with TRISE_NS=200, then 2000, then 5000. Compare
 *   pct_case23_x10, headroom_min, vs_pct_x10_max. If pct_case23 stays ~0 in all runs
 *   the parameter simply never gets used in your operating range.
 */

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "main.h"
#include "mc_type.h"
#include "mc_config.h"
#include "mc_interface.h"
#include "pqd_motor_power_measurement.h"
#include "pwm_curr_fdbk.h"
#include "r3_2_g4xx_pwm_curr_fdbk.h"

/* ---------- configuration ---------- */
#define MYDIAG_MOTOR_IDX     0u           /* M1 */
#define MYDIAG_TIM_CLK_MHZ   170u         /* ADV_TIM_CLK_MHz of the project */
#define MYDIAG_NEG_RPM       100          /* speed below -100 rpm counts as reverse rotation */
#define MYDIAG_ANG_RPM_MIN   2500         /* angle estimate only above this speed (EMF dominates) */
#define MYDIAG_ANG_IQ_MAX    1400         /* ... and below this |Iq| in s16 digits (about 15 A, scale not verified) */
#define MYDIAG_ANG_STEADY_RPM 150         /* window only when |speed - speed_ref| < this (steady state) */
#define MYDIAG_MIN_SAMPLES   200u         /* shorter runs are ignored by MyDiag_CycleEnd */
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
  volatile int32_t  sum_iq;
  volatile int32_t  sum_id;
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
  int16_t  iq_std_s16;       /* ripple (standard deviation) of Iq / Id in the window, digits */
  int16_t  id_std_s16;
  int16_t  iq_mean_s16;      /* mean Iq / Id in the same window (digits) */
  int16_t  id_mean_s16;
  int16_t  pwr_mean_w_x10;   /* mean electrical power [W]x10 in the same window (compare between phase shifts) */
} MyDiagCycle_t;

#define MYDIAG_HIST 8u
volatile MyDiagCycle_t g_myDiagHist[MYDIAG_HIST];
volatile uint32_t g_myDiagHistCount;   /* total number of finished cycles */

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
  g_myDiag.sum_iq2 = 0.0f;
  g_myDiag.sum_id2 = 0.0f;
  g_myDiag.sum_id = 0;

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
      g_myDiag.sum_iq2 += (float)((int32_t)Iqd.q * (int32_t)Iqd.q);
      g_myDiag.sum_id2 += (float)((int32_t)Iqd.d * (int32_t)Iqd.d);
      g_myDiag.sum_id += (int32_t)Iqd.d;
      g_myDiag.sum_pw += PQD_GetAvrgElMotorPowerW(pMPM[MYDIAG_MOTOR_IDX]);
    }
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
    c->id_mean_s16 = (int16_t)(g_myDiag.sum_id / (int32_t)g_myDiag.n_ang);
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
    c->id_mean_s16 = 0;
    c->iq_std_s16 = 0;
    c->id_std_s16 = 0;
  }
  c->n_ang_ms = (g_myDiag.n_ang > 65535u) ? 65535u : (uint16_t)g_myDiag.n_ang;
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
