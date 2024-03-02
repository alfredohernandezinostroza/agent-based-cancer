import sys
import Batch
import DataGenerator
import GraphGenerator
import videoGenerator
import warnings

if __name__ == "__main__":
    if len(sys.argv) == 1:
        raise Exception("Incorrent amount of or unrecognized arguments!")
    if len(sys.argv) >= 3 and sys.argv[1] == "postprocess":
        if sys.argv[2] == "graphic":
            raise Exception("\nError! No known command for postprocessing named \'graphic\'. Did you mean \'graphics\'?")
        if sys.argv[2] == "video":
            raise Exception("\nError! No known command for postprocessing named\'video\'. Did you mean \'videos\'?")
        elif sys.argv[2] != "data":
            raise Exception(f"\nError! No known command for postprocessing named \'{sys.argv[2]}\'. Use \'data\', \'videos\' or \'graphics\'.")
    if len(sys.argv) == 2 or len(sys.argv) == 3:
        raise Exception("Incorrent amount of or unrecognized arguments!")
    elif len(sys.argv) == 4:
        #CORRECT ONES
        if sys.argv[1] == "run":
            total_steps     = int(sys.argv[2])
            interval_steps  = int(sys.argv[3])
            Batch.main_Batch(total_steps, interval_steps)
        elif sys.argv[1] == "postprocess" and sys.argv[2] == "data":
                simulation_folder = sys.argv[3]
                DataGenerator.generate_data(simulation_folder)
        elif sys.argv[1] == "postprocess" and sys.argv[2] == "graphics":
            simulation_folder = sys.argv[3]
            amount_of_pictures = 0 #this will generate all pictures!
            warnings.warn("No amount of pictures given to graphics. Generating all pictures.")
            GraphGenerator.generate_graphs(simulation_folder, amount_of_pictures)
        #ERRORS
        else:
            raise Exception("Incorrent amount of or unrecognized arguments!")
    elif len(sys.argv) == 5:
        if sys.argv[1] == "load":
            simulation_folder = sys.argv[2]
            total_steps     = int(sys.argv[3])
            interval_steps  = int(sys.argv[4])
            Batch.main_Batch(total_steps, interval_steps, simulation_folder)
        elif sys.argv[1] == "postprocess" and sys.argv[2] == "graphics":
            simulation_folder = sys.argv[3]
            amount_of_pictures = sys.argv[4]
            GraphGenerator.generate_graphs(simulation_folder, amount_of_pictures)
        elif sys.argv[1] == "postprocess" and sys.argv[2] == "videos":
            simulation_folder = sys.argv[3]
            frame_rate = sys.argv[4]
            videoGenerator.generate_videos(simulation_folder, frame_rate)
        else:
            raise Exception("Incorrent amount or unrecognized arguments!")