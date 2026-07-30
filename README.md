# Stimulated-echo T2 LUT correction

This toolbox creates a LUT table based on sequence parameters and slice profile using EPG simulations for a range of plausible T2 values in the human brain at 3T. 

This toolbox was developped to correct for T2 maps reconstructed with the GRAPPATINI Siemens sequence: the raw signal is not available and the T2 map is contaminated by stimulated echo times that cannot be corrected with a standard StimFit approach ([Lebel and Wilman, MRM 2010](doi.org/10.1002/mrm.22487)). 

Note on assumptions: 
- This approach assumes a negligible effect from B1         inhomogeneity. Further integration may be necessary to optimize the model for regions subject to strong B1 inhomogeneities. With additional measurements required. 
- The range of T2 values simulated for this LUT table is 30 to 130ms. This range may not match all field strength and tissue and may require adjustments.
- The "GRAPPATINI" mode assumes that the monoexponential fit applied during reconstruction skipped the first echo of the multi echo train length and was reconstructed from complex data, hence the simulation do not account for a noise floor. Applying the model to other sequence schemes may require further adjustments.
- The slice profile needs to be provided or simulated as per [StimFit](https://github.com/usc-mrel/StimFit/tree/main).

## Dependencies
- [StimFit](https://github.com/usc-mrel/StimFit/tree/main)