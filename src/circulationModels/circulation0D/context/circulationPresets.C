/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

Description
    The sub-circuit modules that ship with the library, written in the same
    dictionary syntax a case uses for its own "modules" (see circulationModule).
    A case's module of the same name replaces a preset.

Author
    David Kelly, UCD.

\*---------------------------------------------------------------------------*/

#include "circulationContext.H"
#include "IStringStream.H"

namespace
{

const char* presetText = R"PRESETS(

// ---- Windkessels (unit-agnostic: give the parameters in the case's
// units) ----
// All end on a ground inside the module, so the only terminal is "in".

// C in parallel with R
windkessel2
{
    required    ( R C );
    defaults    { initialPressure 0; }
    components
    {
        C       { type compliance; C @C; initialPressure @initialPressure; }
        R       { type resistor;   R @R; }
        ground  { type ground; }
    }
    connections ( (C R ground) );
    terminals   { in C.in; }
}

// characteristic impedance Rch in front of the two-element
windkessel3
{
    required    ( R Rch C );
    defaults    { initialPressure 0; }
    components
    {
        Rch     { type resistor;   R @Rch; }
        C       { type compliance; C @C; initialPressure @initialPressure; }
        R       { type resistor;   R @R; }
        ground  { type ground; }
    }
    connections ( (Rch C R ground) );
    terminals   { in Rch.in; }
}

// Rch and blood inertance L in series in front of the two-element
windkessel4
{
    required    ( R Rch L C );
    defaults    { initialPressure 0; }
    components
    {
        Rch     { type resistor;   R @Rch; }
        L       { type inductor;   L @L; }
        C       { type compliance; C @C; initialPressure @initialPressure; }
        R       { type resistor;   R @R; }
        ground  { type ground; }
    }
    connections ( (Rch L C R ground) );
    terminals   { in Rch.in; }
}

// ---- Closed-loop 4-chamber human circulation
// (Healy / Korakianitis & Shi) ----
// Written in mmHg, mL, s: the case must use  units mmHg_mL;  Any parameter may
// be overridden, e.g.  heart { type healthyHuman; Emax_lv 2.0; }
healthyHuman
{
    units       mmHg_mL;
    defaults
    {
        V0_lv 5;  P0_lv 1;  Emin_lv 0.1;  Emax_lv 2.5;
        Tes_lv 0.3;  Ted_lv 0.45;  Eshift_lv 0;  Vinit_lv 500;
        V0_rv 10; P0_rv 1;  Emin_rv 0.1;  Emax_rv 1.15;
        Tes_rv 0.3;  Ted_rv 0.45;  Eshift_rv 0;  Vinit_rv 500;
        V0_la 4;  P0_la 1;  Emin_la 0.15; Emax_la 0.25;
        Tes_la 0.045; Ted_la 0.09; Eshift_la 0.92; Vinit_la 20;
        V0_ra 4;  P0_ra 1;  Emin_ra 0.15; Emax_ra 0.25;
        Tes_ra 0.045; Ted_ra 0.09; Eshift_ra 0.92; Vinit_ra 20;

        CQ_av 350; CQ_mv 400; CQ_tv 400; CQ_pv 350;

        Csas 0.08; Rsas 0.003; Lsas 6.2e-5; Psas0 100;
        Csat 1.6;  Rsat 0.05;  Lsat 0.0017; Psat0 100;
        Rsar 0.5;  Rscp 0.52;
        Csvn 20.5; Rsvn 0.075; Psvn0 0;
        Cpas 0.18; Rpas 0.002; Lpas 5.2e-5; Ppas0 30;
        Cpat 3.8;  Rpat 0.01;  Lpat 0.0017; Ppat0 30;
        Rpar 0.05; Rpcp 0.25;
        Cpvn 20.5; Rpvn 0.006; Ppvn0 0;
    }
    components
    {
        LV
        {
            type elastanceChamber; V0 @V0_lv; P0 @P0_lv;
            Emin @Emin_lv; Emax @Emax_lv; Tes @Tes_lv; Ted @Ted_lv;
            Eshift @Eshift_lv; initialVolume @Vinit_lv;
        }
        RV
        {
            type elastanceChamber; V0 @V0_rv; P0 @P0_rv;
            Emin @Emin_rv; Emax @Emax_rv; Tes @Tes_rv; Ted @Ted_rv;
            Eshift @Eshift_rv; initialVolume @Vinit_rv;
        }
        LA
        {
            type elastanceChamber; V0 @V0_la; P0 @P0_la;
            Emin @Emin_la; Emax @Emax_la; Tes @Tes_la; Ted @Ted_la;
            Eshift @Eshift_la; initialVolume @Vinit_la;
        }
        RA
        {
            type elastanceChamber; V0 @V0_ra; P0 @P0_ra;
            Emin @Emin_ra; Emax @Emax_ra; Tes @Tes_ra; Ted @Ted_ra;
            Eshift @Eshift_ra; initialVolume @Vinit_ra;
        }

        AV  { type orificeValve; CQ @CQ_av; }
        MV  { type orificeValve; CQ @CQ_mv; }
        TV  { type orificeValve; CQ @CQ_tv; }
        PV  { type orificeValve; CQ @CQ_pv; }

        SAS { type CRL; C @Csas; R @Rsas; L @Lsas; initialPressure @Psas0; }
        SAT { type CRL; C @Csat; R @Rsat; L @Lsat; initialPressure @Psat0; }
        SAR { type resistor; R @Rsar; }
        SCP { type resistor; R @Rscp; }
        SVN { type CR;  C @Csvn; R @Rsvn; initialPressure @Psvn0; }

        PAS { type CRL; C @Cpas; R @Rpas; L @Lpas; initialPressure @Ppas0; }
        PAT { type CRL; C @Cpat; R @Rpat; L @Lpat; initialPressure @Ppat0; }
        PAR { type resistor; R @Rpar; }
        PCP { type resistor; R @Rpcp; }
        PVN { type CR;  C @Cpvn; R @Rpvn; initialPressure @Ppvn0; }
    }
    connections
    (
        (LV AV SAS SAT SAR SCP SVN RA TV RV PV PAS PAT PAR PCP PVN LA MV LV)
    );
}

// ---- Heart and arteries of the closed-loop pipe study (SI) ----
// outlet side -> C1 -R1-L1- C2 -R2- C3 -R3- atrium -mitral- ventricle
//                                    -aortic(+R4)- C4 -> inlet side
// "in" is the node of C1 (connect the CFD outlet port to it) and "out" the node
// of C4 (connect it to the CFD inlet port):  (outletPort heart inletPort)
// Parameters keep the names of the original study's lumpedParameters.
pipeLoopHeart
{
    units       SI;
    defaults
    {
        C1 1.500e-10; R1 6.666e6; L1 5.333e5; C2 1.500e-9; R2 6.667e8;
        C3 7.501e-9;  R3 1.200e7; C4 7.501e-10;
        Rhat1 5.333e10; Rhat2 5.333e9; R4 1.200e7;
        Ea 9.800e8; Vau 1.000e-6; Pa0 22.66; Ka 4.840e5;
        Ev1 2.467e9; Ev2 -5.600e12; Vvu 4.000e-6; Pv0 120.0; Kv 6.200e4;
        Tsv 0.40; tOffset_v 0.0; Tsa 0.24; tOffset_a 0.544;
        Po 12720; Q1 0; P1 12720; P2 12720; Pi 12720;
        Va 1.4082870749878583e-05; Vv 7.9368207007450098e-05;
    }
    components
    {
        C1 { type compliance; C @C1; initialPressure @Po; }
        R1 { type resistor;   R @R1; }
        L1 { type inductor;   L @L1; initialFlow @Q1; }
        C2 { type compliance; C @C2; initialPressure @P1; }
        R2 { type resistor;   R @R2; }
        C3 { type compliance; C @C3; initialPressure @P2; }
        R3 { type resistor;   R @R3; }
        atrium
        {
            type mixedChamber; V0 @Vau;
            active { E1 @Ea; } passive { P0 @Pa0; K @Ka; }
            activation { type sineSquared; offset @tOffset_a; duration @Tsa; }
            initialVolume @Va;
        }
        mitral { type orificeValve; Rhat @Rhat1; }
        ventricle
        {
            type mixedChamber; V0 @Vvu;
            active { E1 @Ev1; E2 @Ev2; } passive { P0 @Pv0; K @Kv; }
            activation { type sineSquared; offset @tOffset_v; duration @Tsv; }
            initialVolume @Vv;
        }
        aortic { type orificeValve; Rhat @Rhat2; R @R4; }
        C4 { type compliance; C @C4; initialPressure @Pi; }
    }
    connections ( (C1 R1 L1 C2 R2 C3 R3 atrium mitral ventricle aortic C4) );
    terminals   { in C1.in; out C4.out; }
}

)PRESETS";

} // End anonymous namespace


Foam::dictionary Foam::circulationPresets()
{
    IStringStream is(presetText);
    return dictionary(is);
}

// ************************************************************************* //
