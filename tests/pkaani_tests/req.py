import requests
import time

start = time.time()

url = "http://pka-ani-env.eba-rpcu3s9j.us-west-1.elasticbeanstalk.com/upload"
files = {"file": open("mini_hu.pdb", "rb")}
response = requests.post(url, files=files)

print(response.text)

rows = []
pka = eval(response.text)
for key in pka:
    print(key)
    print(pka[key])
    row_dict = dict()
    row_dict["res_num"] = key[1]
    row_dict["res_name"] = pka[key][0]
    row_dict["chain_id"] = key[0]
    row_dict["pKa"] = pka[key][1]
    rows.append(row_dict)


end = time.time()
print(f"{round(end - start, 4)} seconds")
