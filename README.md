# Martian-Lithology

1. Overview
The codebase was developed in a stepwise manner rather than as a single monolithic program. The initial stage focused on testing individual utility and supporting functions. Once these components were established, specialised rock-physics and cementation routines were executed for the different model scenarios. The corresponding inversion scripts were then used for the individual cases, followed by refined plotting and result-generation scripts.
2. Workflow at a Glance
The complete computational workflow can be represented as:
01 BASIC FUNCTIONS
Data handling • density calculations • ensemble utilities • plotting support
↓
02 MODEL BUILDING
Berryman/SCM • cementation • pore/fluid modelling
↓
03 CASE-WISE INVERSION
4 model combinations × 2 scenarios (SC1, SC2)
→ FINAL / REFINED SCRIPTS → PLOTS → RESULTS → INTERPRETATION
3. Stage I — Basic and Supporting Codes
The following files were first run and checked individually. These scripts provide the supporting functions required by the later modelling and inversion stages.
File	Role
assigndatarun	Data assignment and preparation for model runs.
brewermap	Plotting/colour-map support.
bulkdensity	Bulk-density calculation.
bulkdensityweighted	Weighted bulk-density calculation.
crameri	Plotting/colour-map support.
findPartner	Supporting ensemble/partner-search routine.
GetEnsemble	Generation/retrieval of model ensembles.
invCumG	Inversion-related supporting calculation.
labelTrianglePlot	Labelling and formatting of triangle plots.
MoveEnsemble	Ensemble movement/update operation.
myHammer	Supporting model/calculation routine.
sampleG	Generation/handling of G-related samples.
TrianglePlot	Triangle/corner-style parameter visualisation.
4. Stage II — Case-Specific Model Building
After the basic functions were tested, the workflow moved to specialised model components. The repository follows a consistent naming structure in which the model domain, mineral/cement condition, and scenario are encoded in the filename.
4.1 Model Components
myberri — Berryman-related rock-physics calculation.
berryscm — Berryman self-consistent modelling component.
cem — Cementation-related modelling component.
mylogpi — Pore/fluid or log-pi-related modelling component.
4.2 Eight Model / Scenario Combinations
The specialised workflow was applied to four principal model combinations, each evaluated for two scenario variants (SC1 and SC2):
Model domain	Condition	Scenarios
Most-upper	Calcite	SC1, SC2
Most-upper	Halite	SC1, SC2
Upper	Calcite	SC1, SC2
Upper	Halite	SC1, SC2
Thus, the specialised code family is repeated systematically for all eight combinations. For each selected case, the relevant Berryman/SCM, cementation, and pore/fluid routines are used before the corresponding inversion is performed.
5. Stage III — Case-Specific Inversion
Once the physical/model components were prepared for a case, the corresponding main inversion script was executed. The inversion files follow the same naming convention as the model components.
main_inversion_mostupper_calcite_sc1
main_inversion_mostupper_calcite_sc2
main_inversion_mostupper_halite_sc1
main_inversion_mostupper_halite_sc2
main_inversion_upper_calcite_sc1
main_inversion_upper_calcite_sc2
main_inversion_upper_halite_sc1
main_inversion_upper_halite_sc2
Each inversion script corresponds to one specific combination of model domain, condition, and scenario. Keeping these files separate makes the comparison between cases transparent and reproducible.
6. Stage IV — Refined and Final Analysis Scripts
After the case-wise calculations and inversion runs, a set of refined scripts was used for the final analysis and visualisation stage:
File	Purpose
newmyberri	Refined Berryman-related modelling routine.
newcem	Refined cementation routine.
newlogpi	Refined pore/fluid or log-pi routine.
newmaininversion	Refined main inversion workflow.
newplotscript	Final/refined plotting workflow.
plotresultrun	Generation/organisation of result plots.
7. Naming Convention
upper / mostupper — identifies the model/domain grouping.
calcite / halite — identifies the corresponding condition.
sc1 / sc2 — identifies the two scenario variants.
main_inversion — case-specific inversion workflow.
new... — later/refined versions of selected routines.
8. Recommended Workflow
For reproducibility, a user following this repository should read and execute the code conceptually in the following order:
Test the basic/supporting functions individually.
Select one of the four model combinations and its SC1/SC2 scenario.
Run the corresponding Berryman/SCM and cementation components.
Run the corresponding pore/fluid modelling component.
Execute the matching main_inversion script.
Repeat the workflow for the remaining cases as required.
Use the refined scripts for final plots and result comparison.
9. Repository Purpose
This repository serves as a computational record of the internship work. It preserves the progression from individual MATLAB functions to integrated case-wise modelling and inversion. The organisation of the files is intended to make the workflow understandable, reproducible, and suitable for linking the computational results with the methodology and results sections of the internship report.
10. Quick Reference
Basic codes → supporting calculations and utilities
myberri / berryscm / cem / mylogpi → case-specific model construction
main_inversion_* → inversion for each case and scenario
new* / plotresultrun → refined analysis, plots, and final results

