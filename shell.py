import os
import subprocess

def run_another_script():
    current_directory = os.path.dirname(os.path.abspath(__file__))
    script_path = os.path.join(current_directory, 'coord_skyfield.py')
    
    try:
        subprocess.run(['python', script_path])  
    except Exception as e:
        print("Error running the script:", e)

if __name__ == "__main__":
    run_another_script()