import csv
import keyboard
import os

#global variable

def main():
    data_analysis_options = ["Run all", "Cell site Histogram", "Position Histogram"]
    directory_path = "Simulations"
    folder_names = get_folder_names(directory_path)
    selected_option = 0
    update_menu(folder_names, selected_option)
    selected_sim = menu_loop(folder_names, selected_option)
    update_menu(data_analysis_options, selected_option)
    selected_analisys = menu_loop(data_analysis_options, selected_option)

def data_analysis_menu():
    options = ["Create Histogram"]

def update_menu(options, selected_option):
    os.system('cls')
    for i, option in enumerate(options):
        if i == selected_option:
            print("➤", option)
        else:
            print("  ", option)

def menu_loop(options, selected_option):
    while True:
        event = keyboard.read_event()        
        if event.event_type == keyboard.KEY_DOWN:
            if event.name == "down" and selected_option < len(options) - 1:
                selected_option += 1
                update_menu(options, selected_option)
            elif event.name == "up" and selected_option > 0:
                selected_option -= 1
                update_menu(options, selected_option)
            elif event.name == "enter":
                os.system('cls')
                return options[selected_option]
        elif event.event_type == keyboard.KEY_UP and event.name == "esc":
            os.system('cls')
            exit()

def get_folder_names(directory_path):
    folder_names = []
    for folder_name in os.listdir(directory_path):
        if os.path.isdir(os.path.join(directory_path, folder_name)) and folder_name.startswith("Sim"):
            folder_names.append(folder_name)
    return folder_names


if __name__ == "__main__":
    main()