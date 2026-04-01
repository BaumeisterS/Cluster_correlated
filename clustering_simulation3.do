

cd "c:\arbeit\fachgesellschaften\DGoEV\2026\Impuls_ClusterkorrelierterDaten\Daten\"
cap log close 
log using log2.txt, t replace 
/*=================================================================
  dental_simulation_v2.do

  REVIDIERTE VERSION

  Änderungen:
  1. Systematische Variation: Praxen (10–100),
     Patienten/Praxis (5–1000), Zeitpunkte (2–10)
  2. Grafische + tabellarische Darstellung der Szenarien
  3. Exposition: ausschließlich car_new
  4. 50 Bootstrap-Replikationen (Testlauf)
  5. Alle Modelle: Gamma-GLM, log-link
     → exp(β) = relative Änderung im Outcome (PPD)

  Drei-Ebenen-Clusterung:
    Ebene 3: Zahnarztpraxis
    Ebene 2: Patient
    Ebene 1: Zahn × Zeit

  Stata ≥ 16 erforderlich
=================================================================*/

version 16
clear all
set more off
macro drop _all


* ═══════════════════════════════════════════════════════════════
*  TEIL 1 – SIMULATIONSPROGRAMM
* ═══════════════════════════════════════════════════════════════

capture program drop sim_dental
program define sim_dental
    syntax , PRactices(integer) PATients(integer) ///
             TIMEpoints(integer) SEED(integer)    ///
           [ TEeth(integer 28) SAVEas(string) ]

    clear
    set seed `seed'

    local b_car_ppd = 0.60
    local b_car_cal = 0.80
    local b_dur_ppd = 0.20
    local b_dur_cal = 0.25
    local b_car_tl  = 0.70
    local sd_prac   = 0.40
    local sd_pat    = 0.70
    local sd_tooth  = 0.30
    local N = `practices' * `patients'

    quietly {

    /* --- A. PRAXIS-EBENE --- */
    set obs `practices'
    gen int practice_id = _n
    gen double re_prac_c = rnormal(0, `sd_prac')
    gen double re_prac_p = rnormal(0, `sd_prac'*0.85)
    gen double re_prac_t = rnormal(0, `sd_prac'*0.75)
    gen byte prac_urban  = rbinomial(1, 0.60)
    gen byte prac_prev   = rbinomial(1, 0.50)
    tempfile fp
    save `fp'

    /* --- B. PATIENTEN-EBENE --- */
    clear
    set obs `N'
    gen long patient_id  = _n
    gen int  practice_id = ceil(patient_id / `patients')
    gen double z1 = rnormal()
    gen double z2 = rnormal()
    gen double z3 = rnormal()
    gen double re_pat_c = `sd_pat'      *  z1
    gen double re_pat_p = `sd_pat'*0.90 * (0.500*z1 + 0.866*z2)
    gen double re_pat_t = `sd_pat'*0.70 * (0.300*z1 + 0.520*z2 + 0.800*z3)
    drop z1 z2 z3
    gen double age_bl   = max(18, min(75, rnormal(40, 12)))
    gen byte   female   = rbinomial(1, 0.52)
    gen double _u       = runiform()
    gen byte   smoking  = cond(_u<0.50, 0, cond(_u<0.75, 1, 2))
    drop _u
    gen byte   diabetes = rbinomial(1, 0.10)
    gen byte   ses      = ceil(runiform()*5)
    gen byte   ohyg     = floor(runiform()*4)
    gen double bmi      = max(18, min(45, rnormal(26,5)))
    merge m:1 practice_id using `fp', nogen assert(match)

    /* --- C. PATIENT × ZAHN --- */
    expand `teeth'
    bysort patient_id: gen byte tooth_num = _n
    gen byte ttype  = cond(tooth_num<=8,  1,     ///
                      cond(tooth_num<=12, 2,     ///
                      cond(tooth_num<=20, 3, 4)))
    gen byte jaw    = cond(tooth_num<=14, 1, 2)
    gen byte molar  = (ttype==4)
    gen byte post   = (ttype>=3)
    gen double re_tooth = rnormal(0, `sd_tooth')
    egen long  tooth_id = group(patient_id tooth_num)

    /* --- D. PATIENT × ZAHN × ZEIT --- */
    expand `timepoints'
    bysort patient_id tooth_num: gen byte time = _n - 1
    sort patient_id tooth_num time
    gen double age = age_bl + time

    /* --- E. KARIES --- */
    gen double xb_c = cond(time==0, -20,             ///
        -3.00 + 0.010*(age_bl-40) + 0.080*time      ///
        - 0.100*female                               ///
        + 0.350*(smoking==2) + 0.100*(smoking==1)    ///
        + 0.200*diabetes - 0.100*ses + 0.200*ohyg    ///
        + 0.020*(bmi-26)                             ///
        + 0.350*molar + 0.150*(ttype==3)             ///
        + 0.100*(jaw==2)                             ///
        - 0.150*prac_prev + 0.100*prac_urban         ///
        + re_prac_c + re_pat_c + re_tooth)
    gen double h_c = invlogit(xb_c)
    sort patient_id tooth_num time
    by patient_id tooth_num: gen double logS_c = sum(log(1 - h_c))
    gen double S_c = exp(logS_c)
    by patient_id tooth_num: gen double U_c = runiform() if _n==1
    by patient_id tooth_num: replace U_c = U_c[1]
    gen byte has_car = (S_c <= U_c)
    by patient_id tooth_num: gen byte car_new = ///
        (has_car==1 & (has_car[_n-1]==0 | _n==1))
    by patient_id tooth_num: gen byte car_t0 = time if car_new==1
    by patient_id tooth_num: replace  car_t0 = car_t0[_n-1] ///
        if mi(car_t0) & has_car==1
    gen byte dur_car = cond(has_car, time - car_t0, 0)
    drop xb_c h_c logS_c S_c U_c

    /* --- F. PPD / CAL --- */
    gen double mu_ppd = max(0.5,                     ///
        2.50 + `b_car_ppd'*has_car                   ///
        + `b_dur_ppd'*dur_car                        ///
        + 0.008*(age-40)                             ///
        + 0.300*(smoking==2) + 0.120*(smoking==1)    ///
        + 0.200*diabetes - 0.040*ses + 0.120*ohyg    ///
        + 0.200*molar + 0.030*time                   ///
        + 0.30*re_prac_p + 0.40*re_pat_p             ///
        + 0.20*re_tooth)
    gen double ppd = max(1, min(12, rgamma(15, mu_ppd/15)))

    gen double mu_cal = max(0.1,                     ///
        1.20 + `b_car_cal'*has_car                   ///
        + `b_dur_cal'*dur_car                        ///
        + 0.012*(age-40)                             ///
        + 0.400*(smoking==2) + 0.180*(smoking==1)    ///
        + 0.250*diabetes - 0.050*ses + 0.150*ohyg    ///
        + 0.250*molar + 0.050*time                   ///
        + 0.30*re_prac_p + 0.50*re_pat_p             ///
        + 0.25*re_tooth)
    gen double cal = max(0, min(15, rgamma(10, mu_cal/10)))
    drop mu_ppd mu_cal

    /* --- G. PARODONTITIS --- */
    gen byte perio     = (cal >= 3 & ppd >= 4)
    gen byte perio_sev = 0
    replace  perio_sev = 1 if cal>=3 & ppd>=4
    replace  perio_sev = 2 if cal>=4 & ppd>=5
    replace  perio_sev = 3 if cal>=6 & ppd>=5
    bysort patient_id time: egen byte n_perio = total(perio)
    gen byte perio_pat = (n_perio >= 2)
    sort patient_id tooth_num time
    by patient_id tooth_num: gen byte perio_inc = ///
        (perio==1 & (perio[_n-1]==0 | _n==1))

    /* --- H. ZAHNVERLUST --- */
    gen double xb_tl = -5.50                         ///
        + `b_car_tl'*has_car + 0.250*dur_car         ///
        + 1.000*perio + 0.500*(perio_sev>=2)         ///
        + 0.015*(age-40)                             ///
        + 0.250*(smoking==2) + 0.100*(smoking==1)    ///
        + 0.150*diabetes - 0.040*ses + 0.100*ohyg    ///
        + 0.150*molar                                ///
        + re_prac_t + re_pat_t + 0.30*re_tooth
    gen double h_tl = invlogit(xb_tl)
    sort patient_id tooth_num time
    by patient_id tooth_num: gen double logS_tl = sum(log(1 - h_tl))
    gen double S_tl = exp(logS_tl)
    by patient_id tooth_num: gen double U_tl = runiform() if _n==1
    by patient_id tooth_num: replace U_tl = U_tl[1]
    gen byte lost = (S_tl <= U_tl)
    by patient_id tooth_num: gen byte _tln = ///
        (lost==1 & (lost[_n-1]==0 | _n==1))
    by patient_id tooth_num: gen byte lost_t = time if _tln==1
    by patient_id tooth_num: replace  lost_t = lost_t[_n-1] ///
        if mi(lost_t) & lost==1
    foreach v in ppd cal perio perio_sev {
        replace `v' = . if lost==1
    }
    drop xb_tl h_tl logS_tl S_tl U_tl _tln

    /* --- I. ABGELEITETE VARIABLEN --- */
    gen byte atrisk = (lost==0)
    gen byte smk_f  = (smoking==1)
    gen byte smk_c  = (smoking==2)
    compress

    }   // end quietly

    if "`saveas'" != "" {
        save "`saveas'", replace
    }
end


* ═══════════════════════════════════════════════════════════════
*  TEIL 2 – SZENARIOSPEZIFIKATION & ÜBERSICHT
* ═══════════════════════════════════════════════════════════════

/* --- 2A: Szenario-Tabelle erstellen --- */
clear
input byte scen_id int(n_prac n_pat n_time) str12 vary_dim
 1   10   50    5  "Praxen"
 2   20   50    5  "Praxen"
 3   30   50    5  "Praxen"
 4   50   50    5  "Praxen"
 5  100   50    5  "Praxen"
 6   30    5    5  "Patienten"
 7   30   10    5  "Patienten"
 8   30   20    5  "Patienten"
 9   30  100    5  "Patienten"
10   30  200    5  "Patienten"
11   30  500    5  "Patienten"
12   30 1000    5  "Patienten"
13   30   50    2  "Zeitpunkte"
14   30   50    3  "Zeitpunkte"
15   30   50    8  "Zeitpunkte"
16   30   50   10  "Zeitpunkte"
end

gen long n_total   = n_prac * n_pat * 28 * n_time
gen long n_analyse = n_prac * n_pat * 28 * (n_time - 1)
format n_total n_analyse %12.0fc

// Schwellenwerte für rechenintensive Methoden
gen byte do_mixed = (n_analyse <=  500000)
gen byte do_gee   = (n_analyse <=  800000)
gen byte do_boot  = (n_analyse <=  300000)

/* --- 2B: Tabellarische Darstellung --- */
di _n as result "{hline 85}"
di as result "  SIMULATIONSSZENARIEN – ÜBERSICHT"
di as result "{hline 85}"
di as text "  Szenario 3 (30 Praxen, 50 Pat./Praxis, 5 ZP) = Basisszenario"
di as text "  do_mixed/gee/boot = 0: Methode bei diesem N übersprungen"
di as result "{hline 85}"

list scen_id n_prac n_pat n_time vary_dim        ///
     n_total n_analyse do_mixed do_gee do_boot,   ///
     separator(5) noobs abbreviate(12)

save "scenarios.dta", replace

/* --- 2C: Grafische Darstellung --- */

// Graph 1: N Beobachtungen pro Variationsdimension
graph hbar (asis) n_total if vary_dim=="Praxen",            ///
    over(n_prac, label(labsize(small)))                      ///
    bar(1, color(navy%70))                                   ///
    ytitle("N Beobachtungen") title("Praxen variiert")       ///
    subtitle("(50 Pat./Praxis, 5 ZP)")                       ///
    name(gs1, replace) nodraw

graph hbar (asis) n_total if vary_dim=="Patienten",          ///
    over(n_pat, label(labsize(small)))                        ///
    bar(1, color(cranberry%70))                               ///
    ytitle("N Beobachtungen") title("Patienten variiert")     ///
    subtitle("(30 Praxen, 5 ZP)")                             ///
    name(gs2, replace) nodraw

graph hbar (asis) n_total if vary_dim=="Zeitpunkte",         ///
    over(n_time, label(labsize(small)))                       ///
    bar(1, color(forest_green%70))                            ///
    ytitle("N Beobachtungen") title("Zeitpunkte variiert")    ///
    subtitle("(30 Praxen, 50 Pat./Praxis)")                   ///
    name(gs3, replace) nodraw

graph combine gs1 gs2 gs3, rows(1)                            ///
    title("Datensatzgröße nach Szenario")                      ///
    note("28 Zähne pro Patient")
graph export "fig_szenarien_nobs.png", replace width(1400)

// Graph 2: Szenarioraum (Scatter)
twoway                                                         ///
    (scatter n_pat n_prac if vary_dim=="Praxen",               ///
        msym(O) msize(vlarge) mcolor(navy%60))                 ///
    (scatter n_pat n_prac if vary_dim=="Patienten",            ///
        msym(S) msize(vlarge) mcolor(cranberry%60))            ///
    (scatter n_pat n_prac if vary_dim=="Zeitpunkte",           ///
        msym(T) msize(vlarge) mcolor(forest_green%60)),        ///
    yscale(log) ylabel(5 10 20 50 100 200 500 1000)            ///
    xlabel(10 20 30 50 100)                                    ///
    xtitle("Anzahl Praxen") ytitle("Patienten / Praxis (log)") ///
    title("Szenarioraum")                                      ///
    legend(order(1 "Praxen var." 2 "Patienten var."            ///
        3 "Zeitpunkte var.") rows(1) size(small))              ///
    note("Zeitpunkte: 2–10 (im Symbol nicht dargestellt)")
graph export "fig_szenarien_raum.png", replace width(1000)


* ═══════════════════════════════════════════════════════════════
*  TEIL 3 – HAUPTSCHLEIFE: SIMULATION & ANALYSE
*
*  Für jedes Szenario:
*    1) Daten simulieren
*    2) Gamma-GLM (log-link) mit car_new als Exposition
*    3) 6 Methoden: Naiv, Cluster-Patient, Cluster-Praxis,
*       Mixed (meglm/Fallback), GEE, Cluster-Bootstrap (50 Reps)
*    4) ICC aus Nullmodell
*    5) Ergebnisse speichern
* ═══════════════════════════════════════════════════════════════

use "scenarios.dta", clear
local nscen = _N

// Alle Parameter in Locals einlesen (bevor sim_dental clear aufruft)
forvalues i = 1/`nscen' {
    local np_`i'  = n_prac[`i']
    local nn_`i'  = n_pat[`i']
    local nt_`i'  = n_time[`i']
    local dmx_`i' = do_mixed[`i']
    local dge_`i' = do_gee[`i']
    local dbo_`i' = do_boot[`i']
}

// Postfile für Ergebnisse
tempname res
postfile `res'                                              ///
    byte(scen_id) int(n_prac n_pat n_time)                  ///
    long(n_obs n_exp)                                       ///
    double(icc_prac icc_pat)                                ///
    double(b_naive  se_naive)                               ///
    double(b_clpat  se_clpat)                               ///
    double(b_clprac se_clprac)                              ///
    double(b_mixed  se_mixed)                               ///
    double(b_gee    se_gee)                                 ///
    double(b_boot   se_boot)                                ///
    using "results_raw.dta", replace

// ─── Schleife ───
forvalues i = 1/`nscen' {

    local np = `np_`i''
    local nn = `nn_`i''
    local nt = `nt_`i''

    di _n as result "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    di as result "  SZENARIO `i'/`nscen':  " ///
        "`np' Praxen × `nn' Pat. × `nt' ZP"
    di as result "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

    timer clear 1
    timer on 1

    // ── Simulieren ──
    sim_dental, practices(`np') patients(`nn') ///
        timepoints(`nt') seed(`=20240100+`i'') teeth(28)

    // ── Analysedatensatz ──
    keep if atrisk==1 & time > 0
    local nobs = _N
    qui count if car_new == 1
    local nexp = r(N)
    di as text "  N Analysebeobachtungen = " %10.0fc `nobs'
    di as text "  N car_new = 1           = " %10.0fc `nexp'

    // Initialisierung (alles missing)
    local icc_prac = .
    local icc_pat  = .
    foreach m in naive clpat clprac mixed gee boot {
        local b_`m'  = .
        local se_`m' = .
    }

    // Prüfung: genug Exponierte?
    if `nexp' < 20 {
        di as error "  ⚠ Zu wenige Exponierte – Szenario übersprungen"
        post `res' (`i') (`np') (`nn') (`nt') (`nobs') (`nexp') ///
            (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.)
        timer off 1
        continue
    }

    // Kovariaten-Makro
    local X "car_new c.age smk_f smk_c diabetes c.ses c.ohyg molar"

    // ── ICC (Nullmodell, linear) ──
    di as text "  → ICC ..."
    capture mixed ppd || practice_id: || patient_id:, ///
        mle iterate(50) nolog
    if _rc == 0 {
        local v1 = exp(2*_b[lns1_1_1:_cons])
        local v2 = exp(2*_b[lns2_1_1:_cons])
        local ve = exp(2*_b[lnsig_e:_cons])
        local vt = `v1' + `v2' + `ve'
        local icc_prac = `v1' / `vt'
        local icc_pat  = `v2' / `vt'
        di as text "    ICC(Praxis)  = " %6.4f `icc_prac'
        di as text "    ICC(Patient) = " %6.4f `icc_pat'
    }

    // ── M1: Naiv (Gamma log-link, keine Clusterung) ──
    di as text "  → M1 Naiv ..."
    capture glm ppd `X', family(gamma) link(log) nolog
    if _rc == 0 {
        local b_naive  = _b[car_new]
        local se_naive = _se[car_new]
    }

    // ── M2: Cluster-robust SE auf Patientenebene ──
    di as text "  → M2 Cluster-Patient ..."
    capture glm ppd `X', family(gamma) link(log) ///
        vce(cluster patient_id) nolog
    if _rc == 0 {
        local b_clpat  = _b[car_new]
        local se_clpat = _se[car_new]
    }

    // ── M3: Cluster-robust SE auf Praxenebene ──
    di as text "  → M3 Cluster-Praxis ..."
    capture glm ppd `X', family(gamma) link(log) ///
        vce(cluster practice_id) nolog
    if _rc == 0 {
        local b_clprac  = _b[car_new]
        local se_clprac = _se[car_new]
    }

    // ── M4: Mixed (meglm Gamma log-link, Fallback: mixed log(PPD)) ──
    if `dmx_`i'' {
        di as text "  → M4 Mixed ..."
        capture meglm ppd `X' || practice_id: || patient_id:, ///
            family(gamma) link(log) intmethod(laplace)        ///
            iterate(100) nolog
        if _rc == 0 {
            local b_mixed  = _b[car_new]
            local se_mixed = _se[car_new]
            di as text "    (meglm konvergiert)"
        }
        else {
            di as text "    meglm gescheitert – Fallback: mixed ln(PPD)"
            capture {
                gen double log_ppd = ln(ppd)
                mixed log_ppd `X' || practice_id: || patient_id:, ///
                    mle iterate(50) nolog
                local b_mixed  = _b[car_new]
                local se_mixed = _se[car_new]
                drop log_ppd
            }
        }
    }
    else {
        di as text "  → M4 Mixed übersprungen (N > 500.000)"
    }

    // ── M5: GEE (Gamma log-link, exchangeable, robust SE) ──
    if `dge_`i'' {
        di as text "  → M5 GEE ..."
        capture {
            xtset tooth_id time
            xtgee ppd `X', family(gamma) link(log) ///
                corr(exchangeable) vce(robust) nolog
        }
        if _rc == 0 {
            local b_gee  = _b[car_new]
            local se_gee = _se[car_new]
        }
        capture xtset, clear
    }
    else {
        di as text "  → M5 GEE übersprungen (N > 800.000)"
    }

    // ── M6: Cluster-Bootstrap auf Praxenebene (50 Reps) ──
    if `dbo_`i'' {
        di as text "  → M6 Bootstrap (50 Reps) ..."
        capture xtset, clear
        capture noisily bootstrap _b,                        ///
            cluster(practice_id) idcluster(newprac)           ///
            reps(50) seed(`=99900+`i'') nodots:               ///
            glm ppd `X', family(gamma) link(log)
        if _rc == 0 {
            local b_boot  = _b[car_new]
            local se_boot = _se[car_new]
        }
        cap drop newprac
    }
    else {
        di as text "  → M6 Bootstrap übersprungen (N > 300.000)"
    }

    // ── Post ──
    post `res' (`i') (`np') (`nn') (`nt') (`nobs') (`nexp')  ///
        (`icc_prac') (`icc_pat')                              ///
        (`b_naive')  (`se_naive')                             ///
        (`b_clpat')  (`se_clpat')                             ///
        (`b_clprac') (`se_clprac')                            ///
        (`b_mixed')  (`se_mixed')                             ///
        (`b_gee')    (`se_gee')                               ///
        (`b_boot')   (`se_boot')

    timer off 1
    qui timer list 1
    di as text "  ✓ Fertig  (Dauer: " %5.1f r(t1) " Sekunden)"
}

postclose `res'

di _n as result "═══ Alle Szenarien abgeschlossen ═══"


* ═══════════════════════════════════════════════════════════════
*  TEIL 4 – ERGEBNISSE: TABELLEN & GRAFIKEN
* ═══════════════════════════════════════════════════════════════

use "results_raw.dta", clear
merge 1:1 scen_id using "scenarios.dta", keepusing(vary_dim) nogen

/* --- Abgeleitete Größen --- */

// exp(β) = relative Änderung in PPD
foreach m in naive clpat clprac mixed gee boot {
    gen double exp_b_`m' = exp(b_`m')
    gen double ci_lo_`m' = exp(b_`m' - 1.96*se_`m')
    gen double ci_hi_`m' = exp(b_`m' + 1.96*se_`m')
}

// SE-Verhältnisse relativ zum naiven Modell
foreach m in clpat clprac mixed gee boot {
    gen double ratio_`m' = se_`m' / se_naive
}

// Design-Effekt (aus Cluster-Praxis-SE)
gen double deff = (se_clprac / se_naive)^2

// Formate
format exp_b_* ci_lo_* ci_hi_*     %6.3f
format se_*                        %7.4f
format ratio_* deff                %5.2f
format icc_*                       %6.4f
format n_obs n_exp                 %10.0fc

/* ─────────────────────────────────────────────────────────────
   4A.  TABELLARISCHE DARSTELLUNG
   ───────────────────────────────────────────────────────────── */

di _n as result "{hline 85}"
di as result "  ERGEBNISSE: Gamma-GLM (log-link)"
di as result "  Exposition:  car_new  (neue Karies, zahnspezifisch)"
di as result "  Outcome:     PPD  (Probing Pocket Depth)"
di as result "  exp(β):      relative Änderung in PPD"
di as result "{hline 85}"

// Tabelle 1: exp(β) nach Methode
di _n as text "TABELLE 1: exp(β) – Punktschätzer nach Methode"
di as text "  exp(β) > 1: PPD bei neuer Karies höher"
list scen_id n_prac n_pat n_time n_exp                ///
     exp_b_naive exp_b_clpat exp_b_mixed exp_b_gee,    ///
     separator(5) noobs abbreviate(12)

// Tabelle 2: Standardfehler nach Methode
di _n as text "TABELLE 2: Standardfehler von β(car_new)"
list scen_id n_prac n_pat n_time                       ///
     se_naive se_clpat se_clprac se_mixed se_gee se_boot, ///
     separator(5) noobs abbreviate(12)

// Tabelle 3: SE-Verhältnisse und Design-Effekt
di _n as text "TABELLE 3: SE-Verhältnis (Methode / Naiv) und DEFF"
di as text "  Ratio > 1: Naiver SE ist zu klein"
di as text "  DEFF = (SE_cluster-Praxis / SE_naiv)²"
list scen_id n_prac n_pat n_time                       ///
     ratio_clpat ratio_clprac ratio_mixed               ///
     ratio_gee ratio_boot deff,                         ///
     separator(5) noobs abbreviate(12)

// Tabelle 4: ICC
di _n as text "TABELLE 4: Intraklassenkorrelationen (PPD)"
list scen_id n_prac n_pat n_time icc_prac icc_pat,    ///
     separator(5) noobs abbreviate(12)

// Tabelle 5: 95%-Konfidenzintervalle (Basis)
di _n as text "TABELLE 5: exp(β) mit 95%-KI – Basisszenario (Szenario 3)"
foreach m in naive clpat clprac mixed gee boot {
    local eb  = exp_b_`m'[3]
    local lo  = ci_lo_`m'[3]
    local hi  = ci_hi_`m'[3]
    local sem = se_`m'[3]
    di as text "  `m':" _col(22)                        ///
        "exp(β) = " %5.3f `eb'                         ///
        "  [" %5.3f `lo' ", " %5.3f `hi' "]"           ///
        "  SE = " %6.4f `sem'
}

save "results_final.dta", replace


/* ─────────────────────────────────────────────────────────────
   4B.  GRAFISCHE DARSTELLUNG
   ───────────────────────────────────────────────────────────── */

/* === Graph 1: SE-Ratio vs Anzahl Praxen === */
preserve
keep if vary_dim == "Praxen"
sort n_prac
twoway                                                          ///
    (connected ratio_clpat  n_prac,                             ///
        msym(O) lp(solid)     lc(navy)         mc(navy))        ///
    (connected ratio_clprac n_prac,                             ///
        msym(S) lp(dash)      lc(cranberry)    mc(cranberry))   ///
    (connected ratio_mixed  n_prac,                             ///
        msym(T) lp(shortdash) lc(forest_green) mc(forest_green)) ///
    (connected ratio_gee    n_prac,                             ///
        msym(D) lp(longdash)  lc(orange)       mc(orange))     ///
    (connected ratio_boot   n_prac,                             ///
        msym(X) lp(dash_dot)  lc(purple)       mc(purple)),    ///
    yline(1, lc(gs12) lp(dot))                                  ///
    xlabel(10 20 30 50 100)                                     ///
    xtitle("Anzahl Zahnarztpraxen")                             ///
    ytitle("SE-Verhältnis (Methode / Naiv)")                    ///
    title("A: Variation der Praxenanzahl")                      ///
    subtitle("(50 Pat./Praxis, 5 Zeitpunkte)")                  ///
    legend(order(1 "Cluster-Pat." 2 "Cluster-Praxis"           ///
        3 "Mixed" 4 "GEE" 5 "Bootstrap") rows(1) size(vsmall)) ///
    name(g_r1, replace) nodraw
restore

/* === Graph 2: SE-Ratio vs Patienten/Praxis === */
preserve
keep if vary_dim == "Patienten" | scen_id == 3
sort n_pat
twoway                                                          ///
    (connected ratio_clpat  n_pat,                              ///
        msym(O) lp(solid)     lc(navy)         mc(navy))        ///
    (connected ratio_clprac n_pat,                              ///
        msym(S) lp(dash)      lc(cranberry)    mc(cranberry))   ///
    (connected ratio_mixed  n_pat,                              ///
        msym(T) lp(shortdash) lc(forest_green) mc(forest_green)) ///
    (connected ratio_gee    n_pat,                              ///
        msym(D) lp(longdash)  lc(orange)       mc(orange))     ///
    (connected ratio_boot   n_pat,                              ///
        msym(X) lp(dash_dot)  lc(purple)       mc(purple)),    ///
    yline(1, lc(gs12) lp(dot))                                  ///
    xscale(log) xlabel(5 10 20 50 100 200 500 1000)             ///
    xtitle("Patienten pro Praxis (log-Skala)")                  ///
    ytitle("SE-Verhältnis (Methode / Naiv)")                    ///
    title("B: Variation der Clustergröße")                      ///
    subtitle("(30 Praxen, 5 Zeitpunkte)")                       ///
    legend(order(1 "Cluster-Pat." 2 "Cluster-Praxis"           ///
        3 "Mixed" 4 "GEE" 5 "Bootstrap") rows(1) size(vsmall)) ///
    name(g_r2, replace) nodraw
restore

/* === Graph 3: SE-Ratio vs Zeitpunkte === */
preserve
keep if vary_dim == "Zeitpunkte" | scen_id == 3
sort n_time
twoway                                                          ///
    (connected ratio_clpat  n_time,                             ///
        msym(O) lp(solid)     lc(navy)         mc(navy))        ///
    (connected ratio_clprac n_time,                             ///
        msym(S) lp(dash)      lc(cranberry)    mc(cranberry))   ///
    (connected ratio_mixed  n_time,                             ///
        msym(T) lp(shortdash) lc(forest_green) mc(forest_green)) ///
    (connected ratio_gee    n_time,                             ///
        msym(D) lp(longdash)  lc(orange)       mc(orange))     ///
    (connected ratio_boot   n_time,                             ///
        msym(X) lp(dash_dot)  lc(purple)       mc(purple)),    ///
    yline(1, lc(gs12) lp(dot))                                  ///
    xlabel(2 3 5 8 10)                                          ///
    xtitle("Anzahl Messwiederholungen")                         ///
    ytitle("SE-Verhältnis (Methode / Naiv)")                    ///
    title("C: Variation der Zeitpunkte")                        ///
    subtitle("(30 Praxen, 50 Pat./Praxis)")                     ///
    legend(order(1 "Cluster-Pat." 2 "Cluster-Praxis"           ///
        3 "Mixed" 4 "GEE" 5 "Bootstrap") rows(1) size(vsmall)) ///
    name(g_r3, replace) nodraw
restore

/* === Kombinierter Graph === */
graph combine g_r1 g_r2 g_r3, cols(1) imargin(small)           ///
    title("Design-Effekte auf Standardfehler")                   ///
    subtitle("Gamma-GLM (log-link); Exposition = car_new; Outcome = PPD") ///
    note("SE-Verhältnis > 1: Naives Modell unterschätzt den SE" ///
         "Horizontale Linie bei 1 = keine SE-Inflation"         ///
         "50 Bootstrap-Replikationen (Testlauf)")               ///
    xsize(8) ysize(14)
graph export "fig_se_ratio_combined.png", replace width(1200)

/* === Graph 4: Design-Effekt vs Clustergröße === */
use "results_final.dta", clear
preserve
keep if vary_dim == "Patienten" | scen_id == 3
sort n_pat
twoway                                                          ///
    (connected deff n_pat, msym(O) lc(cranberry) mc(cranberry)  ///
        lw(medthick)),                                          ///
    yline(1, lc(gs12) lp(dot))                                  ///
    xscale(log) xlabel(5 10 20 50 100 200 500 1000)             ///
    xtitle("Patienten pro Praxis (log-Skala)")                  ///
    ytitle("Design-Effekt (DEFF)")                              ///
    title("Design-Effekt vs. Clustergröße")                     ///
    subtitle("DEFF = (SE{sub:cluster} / SE{sub:naiv})²")        ///
    note("DEFF > 1: Effektive Stichprobe kleiner als nominal")
graph export "fig_deff_clustersize.png", replace width(1000)
restore

/* === Graph 5: ICC nach Szenario === */
use "results_final.dta", clear
preserve
keep if vary_dim == "Patienten" | scen_id == 3
sort n_pat
twoway                                                          ///
    (connected icc_prac n_pat, msym(O) lc(navy) mc(navy))      ///
    (connected icc_pat  n_pat, msym(S) lc(cranberry) mc(cranberry)), ///
    xscale(log) xlabel(5 10 20 50 100 200 500 1000)             ///
    xtitle("Patienten pro Praxis (log-Skala)")                  ///
    ytitle("ICC")                                               ///
    title("Intraklassenkorrelation vs. Clustergröße")           ///
    legend(order(1 "ICC(Praxis)" 2 "ICC(Patient)") rows(1))
graph export "fig_icc_clustersize.png", replace width(1000)
restore

/* === Graph 6: exp(β) Forest Plot – Basisszenario === */
use "results_final.dta", clear
preserve
keep if scen_id == 3

// Reshape für Forest Plot
gen obs = 1
reshape long exp_b_ ci_lo_ ci_hi_ se_, i(obs) j(method) string

gen byte meth_num = .
replace  meth_num = 6 if method == "naive"
replace  meth_num = 5 if method == "clpat"
replace  meth_num = 4 if method == "clprac"
replace  meth_num = 3 if method == "mixed"
replace  meth_num = 2 if method == "gee"
replace  meth_num = 1 if method == "boot"

label define methlbl                               ///
    6 "Naiv (keine Clusterung)"                    ///
    5 "Cluster-robust (Patient)"                   ///
    4 "Cluster-robust (Praxis)"                    ///
    3 "Mixed Model (meglm)"                        ///
    2 "GEE (exchangeable)"                         ///
    1 "Cluster-Bootstrap (Praxis)"
label values meth_num methlbl

twoway                                                        ///
    (rcap ci_lo_ ci_hi_ meth_num, horizontal lc(navy))        ///
    (scatter meth_num exp_b_, msym(O) mc(cranberry) msize(large)), ///
    xline(1, lc(gs12) lp(dot))                                ///
    xlabel(, format(%5.3f))                                    ///
    xtitle("exp(β): Relative Änderung in PPD")                 ///
    ytitle("")                                                 ///
    ylabel(1/6, valuelabel labsize(small) angle(0))            ///
    title("Basisszenario: exp(β) mit 95%-KI nach Methode")     ///
    subtitle("30 Praxen, 50 Pat./Praxis, 5 Zeitpunkte")       ///
    note("Gestrichelt: exp(β) = 1 (kein Effekt)")             ///
    legend(off)
graph export "fig_forest_base.png", replace width(1000)
restore


* ═══════════════════════════════════════════════════════════════
*  TEIL 5 – ICC & DESIGN-EFFEKTE: BASISSZENARIO VERTIEFT
* ═══════════════════════════════════════════════════════════════

// Basisszenario neu simulieren für detaillierte Analyse
sim_dental, practices(30) patients(50) timepoints(5) ///
    seed(20240103) teeth(28) saveas("sim_base.dta")

use "sim_base.dta", clear
keep if atrisk==1 & time > 0

di _n as result "{hline 70}"
di as result "  VERTIEFUNG: ICC & DESIGN-EFFEKTE (Basisszenario)"
di as result "{hline 70}"

/* --- 5A: ICC aus Nullmodell --- */
di _n as text "─── ICC für PPD (Nullmodell) ───"
mixed ppd || practice_id: || patient_id:, mle nolog
estat icc

/* --- 5B: Varianzkomponenten & Design-Effekt --- */
local v_prac = exp(2*_b[lns1_1_1:_cons])
local v_pat  = exp(2*_b[lns2_1_1:_cons])
local v_res  = exp(2*_b[lnsig_e:_cons])
local v_tot  = `v_prac' + `v_pat' + `v_res'

local icc_prac = `v_prac' / `v_tot'
local icc_pat  = `v_pat'  / `v_tot'
local icc_both = (`v_prac' + `v_pat') / `v_tot'

// Mittlere Clustergrößen
bysort practice_id: gen long N_prac = _N
qui sum N_prac, meanonly
local m_prac = r(mean)

bysort patient_id: gen long N_pat = _N
qui sum N_pat, meanonly
local m_pat = r(mean)

// Design-Effekte
//   DEFF = 1 + (m̄ - 1) × ICC
local deff_prac = 1 + (`m_prac' - 1) * `icc_prac'
local deff_pat  = 1 + (`m_pat'  - 1) * `icc_pat'
local deff_both = 1 + (`m_prac' - 1) * `icc_both'

di _n as text "─── Varianzkomponenten ───"
di as text "  σ²(Praxis)   = " %8.4f `v_prac'
di as text "  σ²(Patient)  = " %8.4f `v_pat'
di as text "  σ²(Residual) = " %8.4f `v_res'

di _n as text "─── ICC ───"
di as text "  ICC(Praxis)          = " %6.4f `icc_prac'
di as text "  ICC(Patient)         = " %6.4f `icc_pat'
di as text "  ICC(Praxis+Patient)  = " %6.4f `icc_both'

di _n as text "─── Mittlere Clustergrößen ───"
di as text "  m̄(Praxis)  = " %8.0f `m_prac'
di as text "  m̄(Patient) = " %8.0f `m_pat'

di _n as text "─── Design-Effekte ───"
di as text "  DEFF(Praxis)          = " %6.2f `deff_prac'
di as text "  DEFF(Patient)         = " %6.2f `deff_pat'
di as text "  DEFF(Praxis+Patient)  = " %6.2f `deff_both'

di _n as text "─── SE-Inflationsfaktoren ───"
di as text "  √DEFF(Praxis)         = " %6.3f sqrt(`deff_prac')
di as text "  √DEFF(Praxis+Patient) = " %6.3f sqrt(`deff_both')

di _n as text "  → Naive SE um Faktor " %4.2f sqrt(`deff_both') ///
    " zu klein,"
di as text "    wenn Clusterung auf beiden Ebenen ignoriert wird."

/* --- 5C: Vollständige Modellergebnisse (Basisszenario) --- */
di _n as text "─── Detaillierte Modellvergleiche ───"
local X "car_new c.age smk_f smk_c diabetes c.ses c.ohyg molar"

// Alle Modelle schätzen und speichern
glm ppd `X', family(gamma) link(log) nolog
estimates store base_naive

glm ppd `X', family(gamma) link(log) vce(cluster patient_id) nolog
estimates store base_clpat

glm ppd `X', family(gamma) link(log) vce(cluster practice_id) nolog
estimates store base_clprac

capture meglm ppd `X' || practice_id: || patient_id:, ///
    family(gamma) link(log) intmethod(laplace) nolog
if _rc == 0 {
    estimates store base_mixed
}
else {
    gen double log_ppd = ln(ppd)
    mixed log_ppd `X' || practice_id: || patient_id:, mle nolog
    estimates store base_mixed
    drop log_ppd
}

xtset tooth_id time
xtgee ppd `X', family(gamma) link(log) corr(exchangeable) ///
    vce(robust) nolog
estimates store base_gee
xtset, clear

bootstrap _b, cluster(practice_id) idcluster(newprac) ///
    reps(50) seed(99903) nodots:                       ///
    glm ppd `X', family(gamma) link(log)
estimates store base_boot
cap drop newprac

// Vergleichstabelle: β (log-Skala)
di _n as text "β-Koeffizienten (log-Skala):"
estimates table base_naive base_clpat base_clprac            ///
    base_mixed base_gee base_boot,                            ///
    keep(car_new) se(%7.4f) b(%7.4f)                          ///
    stats(N) title("Vergleich: β(car_new)")

// exp(β) für alle Methoden
di _n as text "exp(β) – Relative Änderung in PPD bei neuer Karies:"
foreach m in base_naive base_clpat base_clprac                ///
             base_mixed base_gee base_boot {
    estimates restore `m'
    local eb = exp(_b[car_new])
    local se = _se[car_new]
    local lo = exp(_b[car_new] - 1.96*`se')
    local hi = exp(_b[car_new] + 1.96*`se')
    di as text "  `m':" _col(25)                              ///
       "exp(β) = " %5.3f `eb'                                ///
       "  95%-KI [" %5.3f `lo' ", " %5.3f `hi' "]"           ///
       "  SE(β) = " %6.4f `se'
}


* ═══════════════════════════════════════════════════════════════
*  TEIL 6 – ZUSAMMENFASSUNG & INTERPRETATION
* ═══════════════════════════════════════════════════════════════

di _n as result "{hline 70}"
di as result "  ZUSAMMENFASSUNG & INTERPRETATION"
di as result "{hline 70}"

di _n as text "{bf:Modell:}"
di as text "  Gamma-GLM mit Log-Link"
di as text "  log(E[PPD | X]) = β₀ + β₁·car_new + β₂·age + ..."
di as text "  → exp(β₁) = relative Änderung in E[PPD]"
di as text ""
di as text "{bf:Exposition:}"
di as text "  car_new = 1 genau zum Zeitpunkt der neuen Karies"
di as text "  (danach wieder 0 → erfasst akuten Onset-Effekt)"
di as text ""
di as text "{bf:Interpretation von exp(β):}"
di as text "  exp(β) = 1.05:  PPD ist 5% höher bei Zähnen"
di as text "                  mit neuer Karies (zum Onset-Zeitpunkt)"
di as text "  exp(β) = 1.00:  kein Effekt"
di as text ""
di as text "{bf:Erwartete Beobachtungen:}"
di as text ""
di as text "  1. SE-INFLATION BEI IGNORIEREN DER CLUSTERUNG"
di as text "     Ratio SE(Cluster) / SE(Naiv) > 1 in allen Szenarien"
di as text "     Je größer die Cluster (mehr Patienten/Praxis),"
di as text "     desto stärker die SE-Inflation."
di as text "     → Naive p-Werte und KIs sind anti-konservativ!"
di as text ""
di as text "  2. MIXED VS. GEE (Log-Link)"
di as text "     Log-Link ist kollabierbar:"
di as text "     → β(Mixed) ≈ β(GEE) für car_new"
di as text "     → Hauptunterschied in den SE, nicht in exp(β)"
di as text "     (Bei Logit-Link wäre: OR_konditional > OR_marginal)"
di as text ""
di as text "  3. ANZAHL CLUSTER (Praxen)"
di as text "     Wenige Praxen → Sandwich-SE nach unten verzerrt"
di as text "     → Cluster-Bootstrap liefert konservativere SE"
di as text "     Faustregel: g ≥ 30–50 Cluster für zuverlässige"
di as text "     Sandwich-Schätzer"
di as text ""
di as text "  4. DESIGN-EFFEKT"
di as text "     DEFF = 1 + (m̄ - 1) × ICC"
di as text "     Große Cluster + hoher ICC → großer DEFF"
di as text "     → Effektive Stichprobe = N / DEFF"
di as text ""
di as text "{bf:Wahre Simulationsparameter (DGP):}"
di as text "  • Karies → PPD: +0.60 mm  (identity-Skala)"
di as text "  • Kariesdauer → PPD: +0.20 mm/Jahr"
di as text "  • SD: Praxis=0.40, Patient=0.70, Zahn=0.30"
di as text ""
di as text "  Hinweis: car_new erfasst nur den akuten Onset."
di as text "  Im DGP wirkt has_car (kumulativ), daher ist"
di as text "  exp(β_car_new) nicht direkt mit β_DGP=0.60"
di as text "  vergleichbar."

di _n as result "{hline 70}"
di as result "  DATEIEN"
di as result "{hline 70}"
di as text "  scenarios.dta        – Szenariospezifikation"
di as text "  results_raw.dta      – Rohergebnisse (Koeffizienten, SE)"
di as text "  results_final.dta    – Ergebnisse + exp(β), Ratios, DEFF"
di as text "  sim_base.dta         – Basisszenario-Daten"
di as text "  fig_szenarien_*.png  – Szenario-Übersichtsgrafiken"
di as text "  fig_se_ratio_*.png   – SE-Verhältnis-Grafiken"
di as text "  fig_deff_*.png       – Design-Effekt-Grafik"
di as text "  fig_icc_*.png        – ICC-Grafik"
di as text "  fig_forest_*.png     – Forest Plot (Basisszenario)"
di as result "{hline 70}"

exit
log close