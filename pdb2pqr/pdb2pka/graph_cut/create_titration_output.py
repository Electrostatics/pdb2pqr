# TODO: Use pathlib Path instead of os.path
import os.path

def create_output(path, curves):
    file_name_template = "{}_{}_{}.csv"
    for key, curve in curves.items():
        file_name = file_name_template.format(*key)
        file_path = os.path.join(path, file_name)

        with open(file_path, "w") as file:
            for ph_val, value in curve:
                file.write(str(ph_val) + ', ' + str(value) + '\n')
