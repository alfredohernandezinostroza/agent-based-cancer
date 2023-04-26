from Classes.utils import parent_dir, imagesFolder
from Classes.utils import gridsize_utils as gridsize
import cv2
import os

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

def create_video_from_images(images_list, video_path):
    """
    creates a video
    showing the evolution of the images over time.
    """
    # Initialize video writer
    #frame_width, frame_height = cv2.imread(images_list[0]).shape[:2]
    img = cv2.imread(images_list[0])
    height, width, channels = img.shape
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    out = cv2.VideoWriter(video_path, fourcc, 1, (width, height), isColor=True)

    # Add each image to the video
    for image_path in images_list:
        image = cv2.imread(image_path)
        #cv2.imshow(f'image_path:{image_path}', image)

        # If i want i can use enumerate and a if to make the first frame stay longer
        out.write(image)

    # Release video writer
    out.release()
    print("Generated video: " + video_path)



def create_combined_video(images_by_grid_list, video_path):
    """
    Given a list of lists of (timestep, image_path) tuples, creates a video
    showing the evolution of the images over time, with one row for each grid.
    """
    # Initialize video writer
    img = cv2.imread(images_by_grid_list[0][0][1])
    height, width, channels = img.shape
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    out = cv2.VideoWriter(video_path, fourcc, 1, (width, height), isColor=True)

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
    print("Generated video: " + video_path)


if __name__ == "__main__":

    # CHANGE THIS LINE according to the simulation you want to plot the graphs
    nameOfTheSimulation = "Sim maxSteps-14 stepsize-14 N-388 gridsNumber-3"

    # Path where the folder of images are located
    ImagesFolderPath = os.path.join(parent_dir, nameOfTheSimulation, imagesFolder)
    print(f'Analyzing images at {ImagesFolderPath}')

    #Path to save the videos
    VideosFolderPath = os.path.join(parent_dir, nameOfTheSimulation, "Videos")
    print(f'Saving videos at {VideosFolderPath}')

    # Path of the images (.png)
    EcmImagesPath = os.path.join(ImagesFolderPath, "Ecm evolution")
    Mmp2ImagesPath = os.path.join(ImagesFolderPath, "Mmp2 evolution")
    CellsImagesPath = os.path.join(ImagesFolderPath, "Tumor growth")

    # Create folder for all the videos
    if not os.path.exists(VideosFolderPath):
        os.makedirs(VideosFolderPath)
    # If there the visual analysis is already done, tell the user
    else:
        print("This folder already exist! \n\t Videos that already exist will not be overwritten")

    # Create list of lists. Each list will contain all the path to the Ecm of a given grid
    EcmImagesByGrid = (group_images_by_grid(EcmImagesPath))
    Mmp2ImagesByGrid = (group_images_by_grid(Mmp2ImagesPath))
    CellsImagesByGrid = (group_images_by_grid(CellsImagesPath))


    # Create videos for the Ecm of each grid
    for i, images_list in enumerate(EcmImagesByGrid):
        #print(f'i:{i}')
        steps = tuple(image[0] for image in images_list)
        imagesPath = tuple(image[1] for image in images_list)
        #print(f'steps:{steps}')
        #print(f'images:{imagesPath}')
        output_file = os.path.join(VideosFolderPath, f"Ecm Evolution - Grid{i+1}.mp4")
        create_video_from_images(imagesPath, output_file)
    # Create videos for the Ecm for all grids
    output_file = os.path.join(VideosFolderPath, f"Ecm Evolution - All grids.mp4")
    create_combined_video(EcmImagesByGrid, output_file)


    # Create videos for the Mmp2 of each grid
    for i, images_list in enumerate(Mmp2ImagesByGrid):
        #print(f'i:{i}')
        steps = tuple(image[0] for image in images_list)
        imagesPath = tuple(image[1] for image in images_list)
        #print(f'steps:{steps}')
        #print(f'images:{imagesPath}')
        output_file = os.path.join(VideosFolderPath, f"Mmp2 Evolution - Grid{i+1}.mp4")
        create_video_from_images(imagesPath, output_file)
    # Create videos for the Mmp2 for all grids
    output_file = os.path.join(VideosFolderPath, f"Mmp2 Evolution - All grids.mp4")
    create_combined_video(Mmp2ImagesByGrid, output_file)

    
    # Create videos for the Cells of each grid
    for i, images_list in enumerate(CellsImagesByGrid):
        #print(f'i:{i}')
        steps = tuple(image[0] for image in images_list)
        imagesPath = tuple(image[1] for image in images_list)
        #print(f'steps:{steps}')
        #print(f'images:{imagesPath}')
        output_file = os.path.join(VideosFolderPath, f"Cells Evolution - Grid{i+1}.mp4")
        create_video_from_images(imagesPath, output_file)
    # Create videos for the Mmp2 for all grids
    output_file = os.path.join(VideosFolderPath, f"Cells Evolution - All grids.mp4")
    create_combined_video(CellsImagesByGrid, output_file)

print('Finished!')