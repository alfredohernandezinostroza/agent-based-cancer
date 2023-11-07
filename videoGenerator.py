from Classes.utils import parent_dir
from Classes.utils import gridsize_utils as gridsize
import cv2
import os

# To run this code you must be in the parent folder of agent-based-cancer

def group_images_by_grid(images_folder_path):
    """
    Groups all images (.png) on the images_folder_path by grid and returns a list of lists,
    where each list contains the path to the corresponding grid.
    """
    # Initialize dictionary to store paths of images by grid
    images_by_grid = {}

    # Iterate through all images in the folder
    for image_filename in os.listdir(images_folder_path):
        #print(f'image_filename: {image_filename}')
        # Skip non-png files
        if not image_filename.endswith(".png"):
            continue
        
        # Extract grid and timestep from image filename
        grid = image_filename.split("-")[1]
        timestep = image_filename.split("-")[2].split(" ")[0]
        
        # Create list for the grid if it doesn't exist yet
        if grid not in images_by_grid:
            images_by_grid[grid] = []
        
        # Add path to image to the list for the corresponding grid
        image_path = os.path.join(images_folder_path, image_filename)
        images_by_grid[grid].append((timestep, image_path))
    #print(f'IMAGES_BY_GRID: {images_by_grid}')
    # Sort images in each grid list by timestep
    for grid in images_by_grid:
        images_by_grid[grid].sort(key=lambda x: int(x[0].split("step")[-1]))

    # Convert dictionary to list of lists
    images_by_grid_list = [images_by_grid[grid] for grid in sorted(images_by_grid)]

    return images_by_grid_list
# Have to check in the future that it is always ordered


def create_video_from_images(images_list, video_path, frameRate):
    """
    creates a video
    showing the evolution of the images over time.
    """
    # Initialize video writer
    img = cv2.imread(images_list[0])
    height, width, channels = img.shape
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    fps = frameRate
    out = cv2.VideoWriter(video_path, fourcc, fps, (width, height), isColor=True)

    # Add each image to the video
    for image_path in images_list:
        image = cv2.imread(image_path)
        #cv2.imshow(f'image_path:{image_path}', image)

        # If i want i can use enumerate and a if to make the first frame stay longer
        out.write(image)

    # Release video writer
    out.release()
    print(f"\tGenerated video: {video_path}")



def create_combined_video(images_by_grid_list, video_path, frameRate):
    """
    Given a list of lists of (timestep, image_path) tuples, creates a video
    showing the evolution of the images over time, with one row for each grid.
    """
    # Initialize video writer
    img = cv2.imread(images_by_grid_list[0][0][1])
    height, width, channels = img.shape
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    fps = frameRate
    out = cv2.VideoWriter(video_path, fourcc, fps, (width, height), isColor=True)

    # Add each image to the combined video
    for i in range(len(images_by_grid_list[0])):
        # Create a single frame for this timestep
        frame = None
        for j, grid_images in enumerate(images_by_grid_list):
            _, image_path = grid_images[i]
            image = cv2.imread(image_path)
            image = cv2.resize(image, (width // 2, height // 2))
            if frame is None:
                frame = image
            else:
                frame = cv2.hconcat([frame, image])
        # Add frame to combined video
        for aaa in range(3):
            out.write(frame)

        #cv2.imshow('frame',frame)
        #key = cv2.waitKey(0) & 0xFF
        #if key == ord('q'):
        #    break

    # Release video writer
    out.release()
    print(f"\tGenerated video: {video_path}")

def generate_videos(nameOfTheSimulation, frameRate):
    
    # Path where the folder of images are located
    imagesFolder = "Visual analysis"
    ImagesFolderPath = os.path.join(parent_dir, nameOfTheSimulation, imagesFolder)
    print(f'Analyzing images at: {ImagesFolderPath}')

    #Path to save the videos
    VideosFolderPath = os.path.join(parent_dir, nameOfTheSimulation, "Videos")
    print(f'Creating folder to save videos at: {VideosFolderPath}\n')

    print(f'\tFrame rate is: {frameRate}')

    # Path of the images (.png)
    EcmImagesPath = os.path.join(ImagesFolderPath, "Ecm evolution")
    Mmp2ImagesPath = os.path.join(ImagesFolderPath, "Mmp2 evolution")
    CellsImagesPath = os.path.join(ImagesFolderPath, "Tumor growth")
    CellsNumberImagesPath = os.path.join(ImagesFolderPath, "Cells growth")
    VasculatureImagesPath = os.path.join(ImagesFolderPath, "Vasculature evolution")

    # Create folder for all the videos
    if not os.path.exists(VideosFolderPath):
        os.makedirs(VideosFolderPath)
    # If there the visual analysis is already done, tell the user
    else:
        print(f"(Videos folder already exists)! \n\tVideos that already exist will not be overwritten)")

    # Create list of lists. Each list will contain all the path to the Ecm of a given grid
    EcmImagesByGrid = (group_images_by_grid(EcmImagesPath))
    Mmp2ImagesByGrid = (group_images_by_grid(Mmp2ImagesPath))
    CellsImagesByGrid = (group_images_by_grid(CellsImagesPath))
    CellsNumberImagesByGrid = (group_images_by_grid(CellsNumberImagesPath))



    # Create video for Vasculature evolution
    imagesFileNames = [image_filename for image_filename in os.listdir(VasculatureImagesPath)]
    imagesFileNames = sorted(imagesFileNames, key=lambda x: int(x.split('Vasculature-step')[1].split('.')[0]))
    imagesPath = [os.path.join(VasculatureImagesPath, file_name) for file_name in imagesFileNames]
    output_file = os.path.join(VideosFolderPath, f"Vasculature - all grids.mp4")
    create_video_from_images(imagesPath, output_file, frameRate)

    # Create videos for the Ecm of each grid
    for i, images_list in enumerate(EcmImagesByGrid):
        steps = tuple(image[0] for image in images_list)
        imagesPath = tuple(image[1] for image in images_list)
        output_file = os.path.join(VideosFolderPath, f"Ecm Evolution - Grid{i+1}.mp4")
        create_video_from_images(imagesPath, output_file, frameRate)
    # Create videos for the Ecm for all grids
    output_file = os.path.join(VideosFolderPath, f"Ecm Evolution - All grids.mp4")
    create_combined_video(EcmImagesByGrid, output_file, frameRate)

    # Create videos for the Mmp2 of each grid
    for i, images_list in enumerate(Mmp2ImagesByGrid):
        steps = tuple(image[0] for image in images_list)
        imagesPath = tuple(image[1] for image in images_list)
        output_file = os.path.join(VideosFolderPath, f"Mmp2 Evolution - Grid{i+1}.mp4")
        create_video_from_images(imagesPath, output_file, frameRate)
    # Create videos for the Mmp2 for all grids
    output_file = os.path.join(VideosFolderPath, f"Mmp2 Evolution - All grids.mp4")
    create_combined_video(Mmp2ImagesByGrid, output_file, frameRate)

    # Create videos for the Tumor of each grid
    for i, images_list in enumerate(CellsImagesByGrid):
        steps = tuple(image[0] for image in images_list)
        imagesPath = tuple(image[1] for image in images_list)
        output_file = os.path.join(VideosFolderPath, f"Cells Evolution - Grid{i+1}.mp4")
        create_video_from_images(imagesPath, output_file, frameRate)
    # Create videos for the Tumor for all grids
    output_file = os.path.join(VideosFolderPath, f"Cells Evolution - All grids.mp4")
    create_combined_video(CellsImagesByGrid, output_file, frameRate)

    # Create video for the cell number growth
    for i, images_list in enumerate(CellsNumberImagesByGrid):
        steps = tuple(image[0] for image in images_list)
        imagesPath = tuple(image[1] for image in images_list)
        output_file = os.path.join(VideosFolderPath, f"Cells Number - Grid{i+1}.mp4")
        create_video_from_images(imagesPath, output_file, frameRate)
    # Create videos for the Tumor for all grids
    output_file = os.path.join(VideosFolderPath, f"Cells Number - All grids.mp4")
    create_combined_video(CellsNumberImagesByGrid, output_file, frameRate)

    print(f'\nFinished!')


if __name__ == "__main__":

    # CHANGE THIS LINE according to the simulation you want to plot the graphs
    nameOfTheSimulation = "Sim maxSteps-500 stepsize-10 N-388 gridsNumber-3"
    frameRate = 20

    # This runs all the code to generate the videos
    generate_videos(nameOfTheSimulation, frameRate)


