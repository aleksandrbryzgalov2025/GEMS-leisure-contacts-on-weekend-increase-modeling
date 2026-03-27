# Collection of simulation files and results for the case of increased leisure contacts during the weekend
This is a special repository for the simulation results using the **G**erman **E**pidemic **M**icrosimulation **S**ystem (GEMS) version 0.4.3. You'll find an extensive list of tutorials and examples in the official [GEMS documentation](https://immidd.github.io/GEMS/). 
To reproduce the results, the GEMS version 0.4.3 has to be properly installed.
The collection includes the following:  
* Variation in the number of leisure contacts on the weekend day: 0.0; 2.5; 5.0; 7.5; 10.0; 12.5; 15.0; 17.5.
* Variation of the latent period value (leisure contacts are 17.5): 1; 2; 3; 4; 5; 6; 7.
* Overcritical and undercritical cases: latent period and infectious period are both 5 days; in baseline scenario reproduction number is less than 1, in alternative scenario R_0 > 1.
* Week structure variation combined with latent period variation (1 - 7): baseline (uniform contacts daily); 5 working days/ 2 weekend days; 4 working days/ 3 weekend days; 2 working days/ 2 weekend days.
* Week structure variation (baseline (uniform contacts daily); 5 working days/ 2 weekend days; 4 working days/ 3 weekend days; 2 working days/ 2 weekend days) for 2 special cases: latent period equals 4 days and infectious period equals 6 days; latent period equals 6 days and infectious period equals 3 days.
* Additional decrease of the working/ school contacts during the weekend combined with week structure variation (baseline (uniform contacts daily); 5 working days/ 2 weekend days; 4 working days/ 3 weekend days; 2 working days/ 2 weekend days) for 2 special cases: latent period equals 1 day and infectious period equals 5 days; latent period equals 6 days and infectious period equals 3 days.

All results can be found in the [applications folder](applications/). Each subfolder (usually named by lp_ ip_ something, describing a certain simulation) contains 4 files:
* paper_run_...jl - it's a source code file, which should be used to run
* synth_pop_...toml - contains the main parameters of the simulation
* infections_ag.jld2 - contains the aggregated simplified results for 100 runs
* report.pdf - contains aggregated simplified report and infection curves plot

It's not planned to update the code or results, because the present work is finished. 



