import glob
from paraview.simple import *
import os
results_path = "/Users/daby/LargeFiles/rMartin_Results/shape_registration/"
files_RL     = glob.glob(results_path + "Mapping_Fine_mesh_morphed_sphere_PA*_RL_000.vtu")
files_LL     = glob.glob(results_path + "Mapping_Fine_mesh_morphed_sphere_PA*_LL_000.vtu")

files_RL.sort()
files_LL.sort()


path_outputs = "./mappings_lib"

#for mapping in files_RL:
for i in range(len(files_RL)):
   mapping_RL = files_RL[i]
   mapping_LL = files_LL[i]
   reader = XMLUnstructuredGridReader(FileName=[mapping_RL]) 
   reader_LL = XMLUnstructuredGridReader(FileName=[mapping_LL]) 
   renderView = GetActiveViewOrCreate('RenderView')
   renderView.CameraViewUp = [0, 0, -1]
   renderView.CameraFocalPoint = [0, 0, 0]
   renderView.CameraPosition = [0.25, 1, -0.25]
   display = Show(reader, renderView)
   display_LL = Show(reader_LL, renderView)
   warped = WarpByVector(Input=reader)
   warped_LL = WarpByVector(Input=reader_LL)

   warpByVector1Display = Show(warped, renderView, 'UnstructuredGridRepresentation')
   ColorBy(warpByVector1Display, ('POINTS', 'displacement', 'Magnitude'))
   Hide(reader)
   Hide(reader_LL)
   warpByVector1Display = Show(warped_LL, renderView, 'UnstructuredGridRepresentation')
   ColorBy(warpByVector1Display, ('POINTS', 'displacement', 'Magnitude'))
   renderView.ResetCamera()
   renderView.Update()
   
   # get color transfer function/color map for 'displacement'
   displacementLUT = GetColorTransferFunction('displacement')
   
   # Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
   displacementLUT.ApplyPreset('Viridis (matplotlib)', True)   

   # Hide orientation axes
   renderView.OrientationAxesVisibility = 0

   basename = os.path.basename(mapping_RL).replace(".vtu", ".png")
   basename_LL = os.path.basename(mapping_LL).replace(".vtu", ".png")

   SaveScreenshot(path_outputs+basename, renderView)

   Delete(warped)
   Delete(reader)
   Delete(reader_LL)
   Delete(warped_LL)
