/*
 * Configure simulation parameters
 * 
 *
 * Author        : Jorik Stoop
 * Last update   : Nov-2026
 */

#ifndef CONFIG_H
#define CONFIG_H
#include <iomanip>
#include <sstream>
#include <string>

using namespace std;

// Select one model, either wild-type (WT) or a mutated HERG
// #define WT       // wild-type
// #define T474I	// mutated HERG, LQT2
#define R56Q     // mutated HERG, LQT2

// Choose optional drug treatments
// #define ISO_10nM         // isoproterenol at 10nM
// #define K_supplement     // spironolactone
// #define Verapamil_300nM  // verapamil (Calcium channel blocker)

// Choose whether to alter parameters to enable EADs in 1D cable
#define EAD

// Define cycle length as global variable
inline constexpr double BCL = 750.0;    // ms

// Define a tag to add to the output files
inline string result_tag()
{
    string model;
#ifdef WT
    model = "WT";
#endif
#ifdef T474I
    model = "T474I";
#endif
#ifdef R56Q
    model = "R56Q";
#endif

    ostringstream cl;
    cl << fixed << setprecision(2) << BCL;

    string treatment;
#ifdef ISO_10nM
    treatment += "_ISO";
#endif
#ifdef K_supplement
    treatment += "_KSUPP";
#endif
#ifdef Verapamil_300nM
    treatment += "_VERA";
#endif

    if (treatment.empty()){
        treatment = "_NODRUG";
    }
    return model + "_CL_" + cl.str() + treatment;
}

#endif