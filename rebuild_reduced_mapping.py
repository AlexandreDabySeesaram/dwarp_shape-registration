import dolfin
import dolfin_warp as dwarp
import os
import glob

def save_rebuilt_mapping(
        mesh                                                                   ,
        lung                        = "RL"                                     , 
        basename                    = "Barycenter_fine_mapping"                , 
        mappings_basename           = "Mapping_fine_sphere"                    , 
        regul_model                 = "ogdenciarletgeymonatneohookean"         ,
         )
