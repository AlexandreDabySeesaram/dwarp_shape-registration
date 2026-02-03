from PIL import Image
import glob
def compile_images_to_grid(image_paths, grid_size, output_path, scale_factor=0.5):
    # Determine the size of each image
    image = Image.open(image_paths[0])
    image_width, image_height = image.size

    # Create a new image with the size of the grid
    grid_image = Image.new('RGB', (grid_size[1] * image_width, grid_size[0] * image_height))

    # Paste each image into the grid
    for index, path in enumerate(image_paths):
        img = Image.open(path)
        row = index // grid_size[1]
        col = index % grid_size[1]
        grid_image.paste(img, (col * image_width, row * image_height))

    # Resize the grid image
    new_width = int(grid_image.width * scale_factor)
    new_height = int(grid_image.height * scale_factor)
    resized_grid_image = grid_image.resize((new_width, new_height), Image.LANCZOS)

    # Save the resulting grid image
    resized_grid_image.save(output_path)

# List of image paths
#image_paths = glob.glob("/Users/daby/Documents/Code/dwarp_shape-registration/Results/png/selection/*")
image_paths = glob.glob("/Users/daby/Documents/Code/dwarp_shape-registration/post/mappings_lib/*")

# image_paths = ['image1.png', 'image2.png', ..., 'image40.png']  # Replace with your image paths

# Grid size (rows, columns)
grid_size = (6, 6)

# Output path for the compiled grid image
output_path = 'compiled_grid.png'

# Scale factor for resizing
scale_factor = 1  # Adjust this value to resize the output image

# Compile images to grid
compile_images_to_grid(image_paths, grid_size, output_path, scale_factor)
