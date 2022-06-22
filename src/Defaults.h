//
//  Defaults.h
//  SwinMD
//
//  Created by Edoardo Savoia on 15/6/2022.
//

#ifndef Defaults_h
#define Defaults_h

/* FILENAMES */

#define DEFAULT_INFILE "config.txt"

#define DEFAULT_result_FNAME "result.out"
#define DEFAULT_velbehav_FNAME "velbehav.out"
#define DEFAULT_pressureResult_FNAME "pressureResult.out"
#define DEFAULT_pressureAcum_FNAME "pressureAcum.out"
#define DEFAULT_Lustig_FNAME "Lustig.out"
#define DEFAULT_ThermoProp_FNAME "ThermoProp.out"
#define DEFAULT_LustigAverages_FNAME "LustigAverages.out"
#define DEFAULT_resultInduction_FNAME "resultInduction.out"
#define DEFAULT_pTensor_FNAME "pTensor.out"

#define DEFAULT_sysData_FNAME "sysData.out"

/* INTEGRATOR CODES */

#define GearIntegrator_TYPE 0
#define VelocityIntegrator_TYPE 1
#define LFIntegrator_TYPE 2
#define NHIntegrator_TYPE 3

/* CONSTRAINTS CODES */
#define NO_CONSTRAINT 0
#define GAUSSIAN_CONSTRAINT 1
#define SHAKE_CONSTRAINT 2

#endif /* Defaults_h */
