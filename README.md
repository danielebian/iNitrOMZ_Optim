# iNitrOMZ
Requires MATLAB 2013 or above.

### Quick intro
    iNitrOMZ comprises of a biogeochemical model embeded in a below-mixed layer 
    1-D advection diffusion model. Forced with an export flux at the base of the 
    mixed-layer, iNitrOMZ resolves a comprehensive set of processes involved in 
    the remineralization of the sinking organic matter, and especially those that 
    involve N-species.

### Running iNitrOMZ
#### Setting the root path
    Let's call the path to iNitrOMZ/ -- $NITROMSPATH 
      (1) Open the model initialization function
            `$NITROMSPATH/iNitrOMZ_v5.4/UserParams/bgc1d_initialize.m`
            for editing and set: `bgc.root='$NITROMSPATH'`;
      (2) open the template runscript  `$NITROMSPATH/iNitrOMZ_v5.4/runscripts/bgc_run.m`
          for editing and set: `bgc1d_root='$NITROMSPATH/'`
#### Run the model
    Run the template script `$NITROMSPATH/iNitrOMZ_v5.4/runscripts/bgc_run.m` in MATLAB
#### Customizing the run
    Change the model defaults by modifying the initialization scripts 
    in $NITROMSPATH/UserParams/ (see section on Code structure for a detailing of the 
    scripts)



### Code structure
 #### iNitrOMZ_v5.4/runscripts/  
      Template scripts to run or optimize the model
        - bgc_run.m -- template running script
        - gam.m -- template script for optimizing model parameters to observations using 
                   MATLAB's genetic algorithm function (ga() requires MATLAB 2018+ and 
                   the global optimization toolbox)
                   
 #### iNitrOMZ_v5.4/UserParams/ 
      Model user-customizable initialization functions 
        - bgc1d_initialize.m -- main initialization script. The user can modify 
                                general model parameters.
        - bgc1d_initboundary.m -- the user can specify/modify boudary conditions
        - bgc1d_initbgc_params.m -- the user can specify/modify biogeochemical 
                                    parameters
        - bgc1d_initIso_params.m -- the user can specify/modify parameters related 
                                    to N isotopes
  
 #### iNitrOMZ_v5.4/bgc1d_src/
      Core model functions 
        - bgc1d_initialize_DepParam.m -- calculates dependant model parameters
        - bgc1d_initIso_Dep_params.m -- calculates dependant model parameters 
                                        related to N isotopes
        - bgc1d_initIso_update_r15n.m -- used to update N isotopic fractions
        - bgc1d_advection_diff_opt.m -- this is the model core. This function 
                                        performs the advection and diffusion of 
                                        model tracers, applies sources and sinks,
                                        and applies restoring. Also handles model 
                                        output archiving.
        - bgc1d_sourcesink.m -- calculates sources and sinks of model tracers
        - bgc1d_restoring_initialize.m -- initializes lateral restoring forcing
        - bgc1d_restoring.m -- calculates lateral restoring of model tracers
        
 #### iNitrOMZ_v5.4/functions 
      Small functions used for various calculations during intergration
        - not listed here ...
      
 #### iNitrOMZ_v5.4/pprocess/ 
      Post-processing functions usefull for analysing the solution
        - bgc1d_postprocess.m -- processes the final archived model solution into 
                                  a user-friendly structure 
        - more not listed here ...
        
 #### iNitrOMZ_v5.4/Data/
      Data used either for model forcing, or for model validation.
       
 #### iNitrOMZ_v5.4/restart/
      Where restart files are stored
       
 #### iNitrOMZ_v5.4/saveOut/
      Where model output is archived.
>>>>>>> b141387c36af3b771abd895e3698ad0bd2970a81
