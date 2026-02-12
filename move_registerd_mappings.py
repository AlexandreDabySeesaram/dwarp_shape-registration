import os
import shutil

#source_dir = '/Users/daby/LargeFiles/rMartin_Results/shape_registration/' 
#targ_dir = '/Users/daby/Documents/Code/dwarp_shape-registration/Mappings/' 

source_dir = '/Users/daby/LargeFiles/rMartin_Results/shape_registration/' 
targ_dir = '/Users/daby/Documents/Code/dwarp_shape-registration/Mappings/' 


# Create targination directory if it doesn't exist
if not os.path.exists(targ_dir):
    os.makedirs(targ_dir)
lungs = ["LL","RL"]
for lung in lungs:
    for i in range(1, 40):
        # We use :01d or just the number because the original uses 1, 2, 3...
        #old_filename = f"Mapping_Fine_mesh_morphed_sphere_PA{i}_{lung}_000.vtu"
        old_filename = f"Mapping_arm_Fine_mesh_morphed_sphere_PA{i}_{lung}_000.vtu"
        
        # We use :02d to ensure it becomes 01, 02, etc.
        #new_filename = f"Mapping_Fine_mesh_morphed_sphere_{lung}_{i:02d}.vtu"
        new_filename = f"Mapping_arm_Fine_mesh_morphed_sphere_{lung}_{i:02d}.vtu"
        
        source_path = os.path.join(source_dir, old_filename)
        targ_path = os.path.join(targ_dir, new_filename)
        
        if os.path.exists(source_path):
            shutil.move(source_path, targ_path)
           print(f"Moved: {old_filename} -> {new_filename}")
        else:
            print(f"Warning: {old_filename} not found.")

print("Moving all files complete!")
