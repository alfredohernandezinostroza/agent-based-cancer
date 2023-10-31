import csv
import keyboard
import os
import GraphGenerator
import VideoGenerator
import time
import Batch

def main_menu():
    menu_options = ["New simulation", "Load simulation", "Postprocessing", "Exit"]
    default_option = 0
    menu_title = "Welcome to MetaSpread: a Cancer Simulation Program!"
    selected_program = menu_loop(menu_options, default_option, menu_title)
    if selected_program == "New simulation":
        simulation_menu()
    elif selected_program == "Load simulation":
        load_simulation_menu()
    elif selected_program == "Postprocessing":
        postprocessing_menu()
    elif selected_program == "Exit":
        exit()

def simulation_menu():
    print("==================================================================")
    print("New simulation")
    print("==================================================================")
    input()
    
    total_steps = input("Select maximum number of steps: ")
    interval_steps = input("Select the size of the intervals for which the information will be collected: ")
    Batch.main_Batch(int(total_steps), int(interval_steps))    

def load_simulation_menu():
    directory_path = "Simulations"
    folder_names = ["Exit"] + get_folder_names(directory_path)
    selected_option = 0
    postprocessing_title = "Load simulation: select simulation"
    selected_simulation = menu_loop(folder_names, selected_option, postprocessing_title)
    if selected_simulation == "Exit":
        exit()
    input()
    total_steps = input("Select maximum number of steps: ")
    interval_steps = input("Select the size of the intervals for which the information will be collected: ")
    Batch.main_Batch(int(total_steps), int(interval_steps))    
    


def postprocessing_menu():
    directory_path = "Simulations"
    folder_names = ["Exit"] + get_folder_names(directory_path)
    selected_option = 0
    postprocessing_title = "Postprocessing: select simulation"
    selected_simulation = menu_loop(folder_names, selected_option, postprocessing_title)
    if selected_simulation == "Exit":
        exit()
    data_analysis_options = ["Exit", "Run all", "Generate graphs", "Generate videos", "Cell site Histogram", "Position Histogram"]
    while True:
        selected_analysis = menu_loop(data_analysis_options, selected_option, f"Selected option: {selected_simulation}")
        if selected_analysis == "Exit":
            exit()
        if selected_analysis == "Generate graphs":
            GraphGenerator.generate_graphs(selected_simulation)
            time.sleep(3)
        if selected_analysis == "Generate videos":
            frameRate = 20
            VideoGenerator.generate_videos(selected_simulation, frameRate)
            time.sleep(3)
        

def menu_loop(options, selected_option, header = ""):
    update_menu(options, selected_option, header)
    while True:
        event = keyboard.read_event()        
        if event.event_type == keyboard.KEY_DOWN:
            if event.name == "down" and selected_option < len(options) - 1:
                selected_option += 1
                update_menu(options, selected_option, header)
            elif event.name == "up" and selected_option > 0:
                selected_option -= 1
                update_menu(options, selected_option, header)
            elif event.name == "enter":
                os.system('cls')
                return options[selected_option]
        elif event.event_type == keyboard.KEY_UP and event.name == "esc":
            os.system('cls')
            exit()

def update_menu(options, selected_option, header = ""):
    os.system('cls')
    if header != "":
        print("==================================================================")
        print(header)
        print("==================================================================")
    for i, option in enumerate(options):
        if i == selected_option:
            print("➤", option)
        else:
            print("  ", option)

def get_folder_names(directory_path):
    folder_names = []
    for folder_name in os.listdir(directory_path):
        if os.path.isdir(os.path.join(directory_path, folder_name)) and folder_name.startswith("Sim"):
            folder_names.append(folder_name)
    print(folder_names)
    folder_names.sort(key = lambda x: os.path.getmtime(f"{directory_path}\{x}"), reverse = True)
    return folder_names


if __name__ == "__main__":
    main_menu()