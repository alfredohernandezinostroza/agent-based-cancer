import sys
import warnings
import metaspread.simrunner
import metaspread.datagenerator as datagenerator
import metaspread.graphgenerator as graphgenerator
import metaspread.videogenerator as videogenerator

if __name__ == "__main__":
    #simple checks for misspellings in the arguments
    if len(sys.argv) >= 3 and sys.argv[1] == "postprocess":
        if sys.argv[2] == "graphic":
            raise Exception("\nError! No known command for postprocessing named \'graphic\'. Did you mean \'graphics\'?")
        if sys.argv[2] == "video":
            raise Exception("\nError! No known command for postprocessing named\'video\'. Did you mean \'videos\'?")
        elif sys.argv[2] != "data" and sys.argv[2] != "videos" and sys.argv[2] != "graphics" and sys.argv[2] != "all":
            error_string = f"\nError! No known command for postprocessing named \'{sys.argv[2]}\'. Use \'data\', \'videos\' or \'graphics\'."
            raise Exception(error_string)
    #actual processing of the arguments
    if len(sys.argv) == 1:
        from metaspread.interactive import *
        main_menu()
    elif len(sys.argv) == 2:
        if sys.argv[1] == "interactive":
            from metaspread.interactive import *
            main_menu()
        else:
            raise Exception("Incorrent amount of or unrecognized arguments!")
    elif len(sys.argv) == 3:
        raise Exception("Incorrent amount of or unrecognized arguments!")
    elif len(sys.argv) == 4:
        #CORRECT ONES
        if sys.argv[1] == "run":
            total_steps     = int(sys.argv[2])
            interval_steps  = int(sys.argv[3])
            simrunner.main(total_steps, interval_steps)
        elif sys.argv[1] == "postprocess" and sys.argv[2] == "data":
                simulation_folder = sys.argv[3]
                datagenerator.generate_data(simulation_folder)
        elif sys.argv[1] == "postprocess" and sys.argv[2] == "graphics":
            simulation_folder = sys.argv[3]
            amount_of_pictures = 0 #this will generate all pictures!
            warnings.warn("No amount of pictures given to graphics. Generating all pictures.")
            graphgenerator.generate_graphs(simulation_folder, amount_of_pictures)
        #ERRORS
        else:
            raise Exception("Incorrent amount of or unrecognized arguments!")
    elif len(sys.argv) == 5:
        if sys.argv[1] == "load":
            simulation_folder = sys.argv[2]
            total_steps     = int(sys.argv[3])
            interval_steps  = int(sys.argv[4])
            simrunner.main(total_steps, interval_steps, simulation_folder)
        elif sys.argv[1] == "postprocess" and sys.argv[2] == "graphics":
            simulation_folder = sys.argv[3]
            amount_of_pictures = int(sys.argv[4])
            graphgenerator.generate_graphs(simulation_folder, amount_of_pictures)
        elif sys.argv[1] == "postprocess" and sys.argv[2] == "videos":
            simulation_folder = sys.argv[3]
            frame_rate = int(sys.argv[4])
            videogenerator.generate_videos(simulation_folder, frame_rate)
    elif len(sys.argv) == 6:
        if sys.argv[1] == "postprocess" and sys.argv[2] == "all":
                simulation_folder  = sys.argv[3]
                amount_of_pictures = int(sys.argv[4])
                frame_rate = int(sys.argv[5])
                datagenerator.generate_data(simulation_folder)
                datagenerator.generate_data(simulation_folder)
                graphgenerator.generate_graphs(simulation_folder, amount_of_pictures)
                videogenerator.generate_videos(simulation_folder, frame_rate)
        else:
            raise Exception("Incorrent amount or unrecognized arguments!")
    else:
        print(sys.argv)
        raise Exception("Incorrent amount of or unrecognized arguments!")