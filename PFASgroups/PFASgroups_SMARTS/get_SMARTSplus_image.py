from tqdm import tqdm
import requests
import os
import json
import time
import numpy as np
def escape_char(char):
    return char.replace('%','%25').replace('&','%26').replace('+','%2b').replace('#','%23').replace(';','%3B')

def get_smarts_image(smarts, name ='name', filetype = 'svg', write = False):
    """
    visualization_mode 	0 	Visualization mode (0 = Complete Visualization, 1 = ID-Mapping, 2 = Element Symbols, 3 = Structure Diagram-Like)
    visualization_of_default_bonds 	0 	Bond display option (0 = Single Bonds, 1 = Single or Aromatic Bonds)
    labels_for_atoms 	false 	Whether to display labels for atoms
    smartstrim_active 	false 	Enable or disable SMARTS trimming
    legend_mode 	0 	Legend mode (0 = No legend, 1 = Dynamic, 2 = Static, 3 = Both)
    smarts_string_into_picture 	true 	Embed the SMARTS string into the visualization image
    file_format 	svg 	Output file format: either "png" or "svg"
    smileslikearomaticity 	false 	Enable SMILES-like aromaticity interpretation
    detectaromaticbonds 	false 	Detect aromatic bonds in the structure
    """
    apikey = os.environ['SMARTSplus_API_KEY']
    params = {"job_id": hex(hash(smarts)),
                "query": {
                    "smarts": smarts,
                    "parameters": {
                        "file_format": filetype,
                        "visualization_mode": 0,
                        "visualization_of_default_bonds": 1,
                        "labels_for_atoms": False,
                        "smartstrim_active": True,
                        "legend_mode": 3,
                        "smarts_string_into_picture": False,
                        "smileslikearomaticity": False,
                        "detectaromaticbonds": False
                    }
                }
            }
    url = f"https://api.smarts.plus/smartsView/"
    response = requests.post(url, json=params, headers={"X-Api-Key": apikey})
    if response.status_code < 400:
        job_id = response.json()['job_id']
        image = None
        k=0
        status = 0
        while k< 3 and status == 0:
            job_response = requests.get(f"{url}/?job_id={job_id}")
            if job_response.status_code < 400:
                status = 1
                try:
                    image = job_response.json()['result']['image']
                except (KeyError, TypeError) as e:
                    status = 0
                    k+=1
                    time.sleep(np.random.poisson(5))
                    continue
            else:
                k+=1
                time.sleep(np.random.poisson(5))
                continue
        if image is None:
            print(f"Error retrieving SMARTS image for {smarts}")
        if write is True and image is not None:
            with open(f"{name}.{filetype}",'w') as f:
                f.write(image)
        return image
    return None
        
if __name__ == "__main__":
    with open('../data/PFAS_groups_smarts.json','r') as f:
        pfg = json.load(f)
    for g in tqdm(pfg):
        smarts1 = g['smarts1']
        smarts2 = g['smarts2']
        name = g['name']
        id = g['id']
        if smarts1 is not None and smarts1!='':
            if smarts2 is not None and smarts2!='':
                name1 = f'{id} {name} 1'
                x = get_smarts_image(smarts2,name=f'{id} {name} 2', write = True)
            else:
                name1 = f'{id} {name}'
            x = get_smarts_image(smarts1, name = f'{name1}', write = True)