/*   
**   SIR_cvode_py.c
**	 
**   Cvode-C for Python Ctypes implementation of BioNetGen model 'SIR'.
**
**   Code Adapted from templates provided by Mathworks and Sundials.
**
**   Requires the CVODE libraries:  sundials_cvode and sundials_nvecserial.
**   https://computation.llnl.gov/casc/sundials/main.html
**
**-----------------------------------------------------------------------------
**   
**   Compilation notes: 
**   
**   include the model in your C++ file with 
**   
**   #include <SIR.h>
**   
**   and compile with 
**   
**   gcc -fPIC -I/path/to/cvode_lib/include/ -c 'SIR'.c
**   gcc 'SIR'.o --shared -fPIC -L/path/to/cvode_lib/lib/ -lsundials_cvode -lsundials_nvecserial -o 'SIR'.so
**
**   note1: if cvode is in your library path, you can omit path specifications.
** 
**-----------------------------------------------------------------------------
**
**   Usage in Python :
**
**   TODO
*/

/* Library headers */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cvode/cvode.h>             /* prototypes for CVODE  */
#include <nvector/nvector_serial.h>  /* serial N_Vector       */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <cvode/cvode_spgmr.h>       /* prototype for CVSpgmr */

/* Problem Dimensions */
#define __N_PARAMETERS__   3
#define __N_EXPRESSIONS__  5
#define __N_OBSERVABLES__  3
#define __N_RATELAWS__     2
#define __N_SPECIES__      3

// return struct
typedef struct result {
    int status;
    int n_observables;
    int n_species;
    int n_tpts;
    int obs_name_len;
    int spcs_name_len;
    double *observables;
    double *species;
    char *obs_names;
    char *spcs_names;
} RESULT;


/* core function declarations */
RESULT *simulate ( int num_tpts, double *timepts, int num_species_init, double *species_init, int num_parameters, double *parameters );
int   check_flag  ( void *flagvalue, char *funcname, int opt );
void  calc_expressions ( N_Vector expressions, double * parameters );
void  calc_observables ( N_Vector observables, N_Vector species, N_Vector expressions );
void  calc_ratelaws    ( N_Vector ratelaws,  N_Vector species, N_Vector expressions, N_Vector observables );
int   calc_species_deriv ( realtype time, N_Vector species, N_Vector Dspecies, void * f_data );

/* user-defined function declarations */


/* user-defined function definitions  */


/* Calculate expressions */
void
calc_expressions ( N_Vector expressions, double * parameters )
{
    NV_Ith_S(expressions,0) = parameters[0];
    NV_Ith_S(expressions,1) = parameters[1];
    NV_Ith_S(expressions,2) = (1.8/NV_Ith_S(expressions,0));
    NV_Ith_S(expressions,3) = parameters[2];
    NV_Ith_S(expressions,4) = (NV_Ith_S(expressions,0)-NV_Ith_S(expressions,1));
   
}

/* Calculate observables */
void
calc_observables ( N_Vector observables, N_Vector species, N_Vector expressions )
{
    NV_Ith_S(observables,0) = NV_Ith_S(species,0);
    NV_Ith_S(observables,1) = NV_Ith_S(species,1);
    NV_Ith_S(observables,2) = NV_Ith_S(species,2);

}

/* Calculate ratelaws */
void
calc_ratelaws ( N_Vector ratelaws, N_Vector species, N_Vector expressions, N_Vector observables )
{  
    NV_Ith_S(ratelaws,0) = NV_Ith_S(expressions,2)*NV_Ith_S(species,0)*NV_Ith_S(species,1);
    NV_Ith_S(ratelaws,1) = NV_Ith_S(expressions,3)*NV_Ith_S(species,1);

}


/* Calculate species derivatives */
int
calc_species_deriv ( realtype time, N_Vector species, N_Vector Dspecies, void * f_data )
{
    int         return_val;
    N_Vector *  temp_data;
    
    N_Vector    expressions;
    N_Vector    observables;
    N_Vector    ratelaws;

    /* cast temp_data */
    temp_data = (N_Vector*)f_data;
     
    /* sget ratelaws Vector */
    expressions = temp_data[0];
    observables = temp_data[1];
    ratelaws    = temp_data[2];
       
    /* calculate observables */
    calc_observables( observables, species, expressions );
    
    /* calculate ratelaws */
    calc_ratelaws( ratelaws, species, expressions, observables );
                        
    /* calculate derivatives */
    NV_Ith_S(Dspecies,0) = -NV_Ith_S(ratelaws,0);
    NV_Ith_S(Dspecies,1) = NV_Ith_S(ratelaws,0) -NV_Ith_S(ratelaws,1);
    NV_Ith_S(Dspecies,2) = NV_Ith_S(ratelaws,1);


    return(0);
}

/*  Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer 
 */
int check_flag(void *flagvalue, char *funcname, int opt)
{
    int *errflag;

    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (opt == 0 && flagvalue == NULL)
    {
        printf( "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n", funcname );    
        return(1);
    }

    /* Check if flag < 0 */
    else if (opt == 1)
    {
        errflag = (int *) flagvalue;
        if (*errflag < 0)
        {
            printf( "\nSUNDIALS_ERROR: %s() failed with flag = %d\n", funcname, *errflag );
            return(1);
        }
    }

    /* Check if function returned NULL pointer - no memory allocated */
    else if (opt == 2 && flagvalue == NULL)
    {
        printf( "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n", funcname );
        return(1);
    }

    return(0);
}

/*
**   ========================
**   main simulate command 
**   ========================
*/
// we have an input array of (a) time points, (b) species init and (c) parameters
// we have an output array of (a) status, (b) species and (c) observables we need to edit
// plhs = inputs, prhs = outputs, nrhs/nrhs are gone
RESULT *simulate( int num_tpts, double *timepoints, int num_species_init, double *species_init, int num_parameters, double *parameters )
{
    /* variables */
    size_t    i;
    size_t    j;

    /* intermediate data vectors */
    N_Vector  expressions;
    N_Vector  observables;
    N_Vector  ratelaws;

    /* array to hold pointers to data vectors */
    N_Vector  temp_data[3];
    
    /* CVODE specific variables */
    realtype  reltol;
    realtype  abstol;
    realtype  time;
    N_Vector  species;
    void *    cvode_mem;
    int       flag;

    /* make sure timepoints has correct dimensions */
    if ( num_tpts <= 1 )
    {  printf("TIMEPOINTS must be a column vector with 2 or more elements.");  }

    /* make sure species_init has correct dimensions */
    if ( num_species_init != __N_SPECIES__ )
    {  printf("SPECIES_INIT must be a row vector with 3 elements.");  } 

    /* make sure params has correct dimensions */
    if ( num_parameters != __N_PARAMETERS__ )
    {  printf("PARAMS must be a column vector with 3 elements.");  }

    // set output result object
    RESULT res_obj;
    RESULT *result;
    result = malloc(sizeof(RESULT));
    res_obj.n_observables = __N_OBSERVABLES__;
    res_obj.n_species = __N_SPECIES__;
    res_obj.n_tpts = num_tpts;
    double *res_species_ptr;
    double *res_observables_ptr;
    res_species_ptr = malloc(num_tpts*__N_SPECIES__*sizeof(double));
    res_observables_ptr = malloc(num_tpts*__N_OBSERVABLES__*sizeof(double));
    res_obj.species = res_species_ptr;
    res_obj.observables = res_observables_ptr;
    res_obj.status = 0;
    res_obj.obs_name_len = 0;
    res_obj.spcs_name_len = 0;
    // store the observable and species names
    char onames[] = "S/I/R/";
    char snames[] = "S()/I()/R()/";
    char *ptr_obs_names = strdup(onames);
    char *ptr_spcs_names = strdup(snames);
    size_t olen = sizeof(onames);
    size_t slen = sizeof(snames);
    res_obj.obs_names = ptr_obs_names;
    res_obj.spcs_names = ptr_spcs_names;
    res_obj.obs_name_len = olen;
    res_obj.spcs_name_len = slen;

    *result = res_obj;

    /* initialize intermediate data vectors */
    expressions  = NULL;
    expressions = N_VNew_Serial(__N_EXPRESSIONS__);
    if (check_flag((void *)expressions, "N_VNew_Serial", 0))
    {
        result->status = 1;
        return result;
    }

    observables = NULL;
    observables = N_VNew_Serial(__N_OBSERVABLES__);
    if (check_flag((void *)observables, "N_VNew_Serial", 0))
    {
        N_VDestroy_Serial(expressions);
        result->status = 1;
        return result;
    }

    ratelaws    = NULL; 
    ratelaws = N_VNew_Serial(__N_RATELAWS__);
    if (check_flag((void *)ratelaws, "N_VNew_Serial", 0))
    {   
        N_VDestroy_Serial(expressions);
        N_VDestroy_Serial(observables);        
        result->status = 1;
        return result;
    }

    /* set up pointers to intermediate data vectors */
    temp_data[0] = expressions;
    temp_data[1] = observables;
    temp_data[2] = ratelaws;

    /* calculate expressions (expressions are constant, so only do this once!) */
    calc_expressions( expressions, parameters );

    /* SOLVE model equations! */
    species   = NULL;
    cvode_mem = NULL;

    /* Set the scalar relative tolerance */
    reltol = 1e-08;
    abstol = 1e-06;

    /* Create serial vector for Species */
    species = N_VNew_Serial(__N_SPECIES__);
    if (check_flag((void *)species, "N_VNew_Serial", 0))
    {  
        N_VDestroy_Serial(expressions);
        N_VDestroy_Serial(observables);
        N_VDestroy_Serial(ratelaws);
        result->status = 1;
        return result;
    }
    for ( i = 0; i < __N_SPECIES__; i++ )
    {   NV_Ith_S(species,i) = species_init[i];   }
    
    /* write initial species populations into species_out */
    for ( i = 0; i < __N_SPECIES__; i++ )
    {   res_species_ptr[i*num_tpts] = species_init[i];   }
    
    /* write initial observables populations into species_out */ 
    calc_observables( observables, species, expressions );  
    for ( i = 0; i < __N_OBSERVABLES__; i++ )
    {   res_observables_ptr[i*num_tpts] = NV_Ith_S(observables,i);   }

    /*   Call CVodeCreate to create the solver memory:    
     *   CV_ADAMS or CV_BDF is the linear multistep method
     *   CV_FUNCTIONAL or CV_NEWTON is the nonlinear solver iteration
     *   A pointer to the integrator problem memory is returned and stored in cvode_mem.
     */
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    if (check_flag((void *)cvode_mem, "CVodeCreate", 0))
    {                                  
        N_VDestroy_Serial(expressions);
        N_VDestroy_Serial(observables);
        N_VDestroy_Serial(ratelaws);   
        N_VDestroy_Serial(species);    
        CVodeFree(&cvode_mem);         
        result->status = 1;            
        return result;                 
    }                                  



    /*   Call CVodeInit to initialize the integrator memory:     
     *   cvode_mem is the pointer to the integrator memory returned by CVodeCreate
     *   rhs_func  is the user's right hand side function in y'=f(t,y)
     *   T0        is the initial time
     *   y         is the initial dependent variable vector
     */
    flag = CVodeInit(cvode_mem, calc_species_deriv, timepoints[0], species);
    if (check_flag(&flag, "CVodeInit", 1))
    {                                  
        N_VDestroy_Serial(expressions);
        N_VDestroy_Serial(observables);
        N_VDestroy_Serial(ratelaws);   
        N_VDestroy_Serial(species);    
        CVodeFree(&cvode_mem);         
        result->status = 1;            
        return result;                 
    }                                  
   
    /* Set scalar relative and absolute tolerances */
    flag = CVodeSStolerances(cvode_mem, reltol, abstol);
    if (check_flag(&flag, "CVodeSStolerances", 1))
    {                                  
        N_VDestroy_Serial(expressions);
        N_VDestroy_Serial(observables);
        N_VDestroy_Serial(ratelaws);   
        N_VDestroy_Serial(species);    
        CVodeFree(&cvode_mem);         
        result->status = 1;            
        return result;                 
    }                                                                       
   
    /* pass params to rhs_func */
    flag = CVodeSetUserData(cvode_mem, &temp_data);
    if (check_flag(&flag, "CVodeSetFdata", 1))
    {                                  
        N_VDestroy_Serial(expressions);
        N_VDestroy_Serial(observables);
        N_VDestroy_Serial(ratelaws);   
        N_VDestroy_Serial(species);    
        CVodeFree(&cvode_mem);         
        result->status = 1;            
        return result;                 
    }                                  
    
    /* select linear solver */
    flag = CVDense(cvode_mem, __N_SPECIES__);
    if (check_flag(&flag, "CVDense", 1))
    {                                  
        N_VDestroy_Serial(expressions);
        N_VDestroy_Serial(observables);
        N_VDestroy_Serial(ratelaws);   
        N_VDestroy_Serial(species);    
        CVodeFree(&cvode_mem);         
        result->status = 1;            
        return result;                 
    }                                                                    
    
    flag = CVodeSetMaxNumSteps(cvode_mem, 2000);
    if (check_flag(&flag, "CVodeSetMaxNumSteps", 1))
    {                                  
        N_VDestroy_Serial(expressions);
        N_VDestroy_Serial(observables);
        N_VDestroy_Serial(ratelaws);   
        N_VDestroy_Serial(species);    
        CVodeFree(&cvode_mem);         
        result->status = 1;            
        return result;                 
    }                                                                    

    flag = CVodeSetMaxErrTestFails(cvode_mem, 7);
    if (check_flag(&flag, "CVodeSetMaxErrTestFails", 1))
    {                                  
        N_VDestroy_Serial(expressions);
        N_VDestroy_Serial(observables);
        N_VDestroy_Serial(ratelaws);   
        N_VDestroy_Serial(species);    
        CVodeFree(&cvode_mem);         
        result->status = 1;            
        return result;                 
    }                                                                    

    flag = CVodeSetMaxConvFails(cvode_mem, 10);
    if (check_flag(&flag, "CVodeSetMaxConvFails", 1))
    {                                  
        N_VDestroy_Serial(expressions);
        N_VDestroy_Serial(observables);
        N_VDestroy_Serial(ratelaws);   
        N_VDestroy_Serial(species);    
        CVodeFree(&cvode_mem);         
        result->status = 1;            
        return result;                 
    }                                                                    

    flag = CVodeSetMaxStep(cvode_mem, 0.0);
    if (check_flag(&flag, "CVodeSetMaxStep", 1))
    {                                  
        N_VDestroy_Serial(expressions);
        N_VDestroy_Serial(observables);
        N_VDestroy_Serial(ratelaws);   
        N_VDestroy_Serial(species);    
        CVodeFree(&cvode_mem);         
        result->status = 1;            
        return result;                 
    }                                                                    

    /* integrate to each timepoint */
    for ( i=1;  i < num_tpts;  i++ )
    {
        flag = CVode(cvode_mem, timepoints[i], species, &time, CV_NORMAL);
        if (check_flag(&flag, "CVode", 1))
        {
            N_VDestroy_Serial(expressions);
            N_VDestroy_Serial(observables);           
            N_VDestroy_Serial(ratelaws);
            N_VDestroy_Serial(species);
            CVodeFree(&cvode_mem);
            result->status = 1; 
            return result;
        }

        /* copy species output from nvector to matlab array */
        for ( j = 0; j < __N_SPECIES__; j++ )
        {   res_species_ptr[j*num_tpts + i] = NV_Ith_S(species,j);   }
        
        /* copy observables output from nvector to matlab array */
        calc_observables( observables, species, expressions );         
        for ( j = 0; j < __N_OBSERVABLES__; j++ )
        {   res_observables_ptr[j*num_tpts + i] = NV_Ith_S(observables,j);   }      
    }
 
    /* Free vectors */
    N_VDestroy_Serial(expressions);
    N_VDestroy_Serial(observables);  
    N_VDestroy_Serial(ratelaws);        
    N_VDestroy_Serial(species);

    /* Free integrator memory */
    CVodeFree(&cvode_mem);

    return result;
}

void free_result(RESULT *r) {
    free(r->obs_names);
    free(r->spcs_names);
    free(r->observables);
    free(r->species);
    free(r);
}
