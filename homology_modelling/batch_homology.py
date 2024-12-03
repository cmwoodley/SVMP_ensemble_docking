import requests
from urllib import request
import os
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-f","--fasta",type = str, help="Fasta format file to build models for.")
parser.add_argument("-t","--token",help="SWISS-MODEL API Token")

args = parser.parse_args()
token = args.token
file = args.fasta

with open(file,"r") as f:
    fasta = f.readlines()
fasta = [line.strip() for line in fasta]

seq_inds = []
for i, line in enumerate(fasta):
    if line.find(">") != -1: seq_inds+=[i]
seq_inds += [len(fasta)]

seq_uniprot = [fasta[j].split("|")[1] for j in seq_inds[:-1]]
seq_titles = [fasta[j].split("_")[-1] for j in seq_inds[:-1]]

full_seq = []
for i,ind in enumerate(seq_inds[:-1]):
    full_seq += ["".join(fasta[seq_inds[i]+1:seq_inds[i+1]]).strip()]

project_ids = []

if os.path.isdir("./models/") != True:
    os.mkdir("./models/")

for i in range(len(full_seq)):

    if full_seq[i].find("X") != -1: 
        print(f"Sequence { i } contains X")
        continue
    
    response = requests.post(
        "https://swissmodel.expasy.org/automodel",
        headers={ "Authorization": f"Token {token}" },
        json={ 
            "target_sequences": 
                [
                    full_seq[i]
                ],
            "project_title":seq_titles[i]
        })
    
    # Obtain the project_id from the response created above
    project_id = response.json()["project_id"]
    project_ids += [project_id]
    # And loop until the project completes
    import time
    while True:
        # We wait for some time
        time.sleep(15)

        # Update the status from the server 
        response = requests.get(
            f"https://swissmodel.expasy.org/project/{ project_id }/models/summary/", 
            headers={ "Authorization": f"Token {token}" })

        # Update the status
        status = response.json()["status"]

        print('Job status is now', status)

        if status in ["COMPLETED", "FAILED"]:
            break

    if status == "COMPLETED":
        path = "./models/{}".format(seq_uniprot[i])
        os.mkdir(path)

        report = [seq_titles[i]+"\n"]
        for model in response.json()["models"]:
            report += ["model_ID: {}\n".format(model["model_id"]), "GMQE: {}\n".format(model["gmqe"]),"\n"]
            url = model["coordinates_url"]
            request.urlretrieve(url,path+"/"+seq_uniprot[i]+"_"+url.split("/")[-1])

        with open(path+"/report.txt", "w") as f:
            for line in report:
                f.writelines(line)
            