from pynput import keyboard
import os
import metaspread.datagenerator
import metaspread.graphgenerator
import metaspread.videogenerator
import time
import metaspread.simrunner as simrunner

selected_option_index = 0
selected_option = ""
options_list =  ["New simulation", "Load simulation", "Postprocessing", "Exit"]
banner_message = "Welcome to MetaSpread: a Cancer Simulation Program!"

def main_menu():
    with keyboard.Listener(on_press = on_press) as listener:
        os.system('cls')
        print_menu()
        listener.join()
        if selected_option == "New simulation":
            simulation_menu()
        elif selected_option == "Load simulation":
            load_simulation_menu()
        elif selected_option == "Postprocessing":
            postprocessing_menu()
        elif selected_option == "Exit":
            os._exit(0)


def simulation_menu():
    print("==================================================================")
    print("New simulation")
    print("==================================================================")
    total_steps = get_int_input("Select the total timesteps the simulation will take: ") #convert this to time and show, asking if it;s ok instead of the other comments
    interval_steps = get_int_input("Select the size of the intervals for which the information will be collected: ")
    simrunner.main(total_steps, interval_steps)


def load_simulation_menu():
    global selected_option, selected_option_index, banner_message, options_list
    directory_path = "Simulations"
    folder_names = ["Exit"] + get_folder_names(directory_path)
    options_list = folder_names
    selected_option_index = 0
    banner_message = "Load simulation: select simulation"
    with keyboard.Listener(on_press = on_press) as listener:
        os.system('cls')
        print_menu()
        listener.join()
        if selected_option == "Exit":
            os._exit(0)
    total_steps = get_int_input("Select additional number of steps until the loaded simulation will continue: ")
    interval_steps = get_int_input("Select the size of the intervals for which the information will be collected: ")
    print("Loading simulation...")
    simrunner.main(total_steps, interval_steps, loadedSimulationPath=os.path.join(os.getcwd(),directory_path, selected_option))

def postprocessing_menu():
    global selected_option, selected_option_index, banner_message, options_list
    directory_path = "Simulations"
    folder_names = ["Exit"] + get_folder_names(directory_path)
    options_list = folder_names
    selected_option_index = 0
    banner_message = "Postprocessing: select simulation"
    with keyboard.Listener(on_press = on_press) as listener:
        os.system('cls')
        print_menu()
        listener.join()
        if selected_option == "Exit":
            os._exit(0)
    options_list = ["Exit", "Run all", "Begin data analysis", "Begin graphical analysis", "Generate videos"]
    selected_simulation = selected_option
    selected_option_index = 0
    banner_message = f"Postprocessing simulation at {selected_simulation}"
    with keyboard.Listener(on_press = on_press) as listener:
        os.system('cls')
        print_menu()
        listener.join()
        if selected_option == "Exit":
            os._exit(0)
        if selected_option == "Run all":
            amount_of_pictures = get_int_input(f"Select the amount of pictures that will be produced for the graphical analysis: ")
            frame_rate = get_int_input("Select the framerate for the video.\nA framerate of 20 is suggested \nfor large simulations: ")
            datagenerator.generate_data(selected_simulation)
            time.sleep(3)
            graphgenerator.generate_graphs(selected_simulation, amount_of_pictures)
            time.sleep(3)
            videogenerator.generate_videos(selected_simulation, frame_rate)
            time.sleep(3)
        if selected_option == "Begin data analysis":
            datagenerator.generate_data(selected_simulation)
            time.sleep(3)
        if selected_option == "Begin graphical analysis":
            graphical_analysis_menu(selected_simulation)
        if selected_option == "Generate videos":
            frame_rate = get_int_input("Select the framerate for the video.\nA framerate of 20 is sugggested \nfor large simulations:")
            videogenerator.generate_videos(selected_simulation, frame_rate)
            time.sleep(3)



def graphical_analysis_menu(selected_simulation):
    global selected_option, selected_option_index, banner_message, options_list
    directory_path = "Simulations"
    options_list = ["Generate 5 picture summary", "Custom amount of pictures", "Generate All"]
    selected_option_index = 0
    banner_message = "Graphical analysis: select which frames are going to be plotted"
    with keyboard.Listener(on_press = on_press) as listener:
        os.system('cls')
        print_menu()
        listener.join()
        if selected_option == "Exit":
            os._exit(0)
        if selected_option == "Generate 5 picture summary":
            graphgenerator.generate_graphs(selected_simulation, 5)
            time.sleep(3)
        if selected_option == "Custom amount of pictures":
            amount_of_pictures = get_int_input(f"Select the amount of pictures that will be produced: ")
            graphgenerator.generate_graphs(selected_simulation, amount_of_pictures)
            time.sleep(3)
        if selected_option == "Generate All":
            graphgenerator.generate_graphs(selected_simulation)
            time.sleep(3)


def print_menu():
    if banner_message != "":
        print("==================================================================")
        print(banner_message)
        print("==================================================================")
    for i, option in enumerate(options_list):
        if i == selected_option_index:
            print("➤", option)
        else:
            print("  ", option)

def on_press(key):
    global selected_option_index
    if key == keyboard.Key.down:
        if selected_option_index < len(options_list) - 1:
            selected_option_index += 1
            os.system('cls')
            print_menu()
    elif key == keyboard.Key.up:
        if selected_option_index > 0:
            selected_option_index -= 1
            os.system('cls')
            print_menu()
    elif key == keyboard.Key.enter:
        global selected_option
        selected_option = options_list[selected_option_index]
        os.system('cls')
        return False
    elif key == keyboard.Key.esc:
        os._exit(0)

def get_folder_names(directory_path):
    folder_names = []
    for folder_name in os.listdir(directory_path):
        if os.path.isdir(os.path.join(directory_path, folder_name)) and (folder_name.startswith("Sim") or folder_name.startswith("Continu")):
            folder_names.append(folder_name)
    print(folder_names)
    folder_names.sort(key=lambda x: os.path.getmtime(f"{directory_path}\{x}"), reverse=True)
    return folder_names

def get_int_input(text):
    result = ''
    while not isinstance(result, int) or result < 0:
        while result == '':
            os.system('cls')
            result = input(text)
            try:
                result = int(result)
            except:
                result = ''
    return result

if __name__ == "__main__":
    main_menu()