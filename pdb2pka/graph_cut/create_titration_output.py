import os.path

def create_output(path, curves):
    file_name_template = "{}_{}_{}.csv"
    for key, curve in curves.items():
        file_name = file_name_template.format(*key)
        file_path = os.path.join(path, file_name)

        with open(file_path, "w") as f:
            for pH, value in curve:
                f.write(str(pH)+', '+str(value)+'\n')