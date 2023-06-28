import requests
import json
from bs4 import BeautifulSoup
from rdkit import Chem

import pandas as pd

def get_smi_json():
    url = "http://www.dcaiku.com/v1/perseus/mol_search"
    # url = "http://www.dcaiku.com/v1/account"
    url_get_token = "http://www.dcaiku.com/v1/public/token_refresh"

    authorization_token = "eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJzdWIiOiJhY2Nlc3NUb2tlbiIsImNyZWF0ZWRBdCI6IjE2ODc5MjkyNTMxODgiLCJleHBpcmVkQXQiOiIxNjg3OTM2NDUzMTg4Iiwicm9sZSI6IlJPTEVfVVNFUiIsImp0aSI6InBlcnNldXMiLCJlbWFpbCI6Ijk3MjMxMzAwNkBxcS5jb20ifQ.U1umFXzhd1abS8458TZ3aAVQY0MUFCxRGRBphykzAjI"
    headers = {
        "Accept": "application/json, text/plain, */*",
        "Accept-Encoding": "gzip, deflate",
        "Accept-Language": "zh-CN,zh;q=0.9",
        "Authorization": authorization_token,
        "Content-Length": "36",
        "Content-Type": "application/json;charset=UTF-8",
        "Host": "www.dcaiku.com",
        "Origin": "http://www.dcaiku.com",
        "Proxy-Connection": "keep-alive",
        "Referer": "http://www.dcaiku.com/",
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/114.0.0.0 Safari/537.36"
    }

    data_json = [
        {"label":"c","type":"title"}
    ]
    
    refresh_json = {
        "refreshToken": authorization_token
    }
    k_count = 0
    with open("data/mol_id_6_28.txt", "w") as f:
        while True:
            try:
                r = requests.post(url, data=json.dumps(data_json), headers=headers)
                data = json.loads(r.text)
                
                for j in range(len(data["data"]["moleculeItemVos"])):
                    f.write(data["data"]["moleculeItemVos"][j]["id"] + "\n")
                
                data_json[0]["label"] = data["data"]["next"]
                data_json[0]["type"] = "next"
                
                r = requests.post(url_get_token, data=json.dumps(refresh_json), headers=headers)
                token_json = json.loads(r.text)
                if "error" in token_json.keys():
                    print(r.text)
                    break
                authorization_token = token_json["message"]
                headers["Authorization"] = authorization_token
                refresh_json["refreshToken"] = authorization_token
                
                k_count += 1
                
            except:
                print(k_count)
                continue
    
# get_smi_json()

def get_graph_matrix(smi):
    url = "http://cimg.dcaiku.com/"
    
    data_json = {
        "csrfmiddlewaretoken": "rJYBYZoM4pLEIG653my8aw8YfJPYvm9SSlsrc8gmZlVhvuKcqdfOo6jIdQibj6sJ",
        # "smiles": "c1cc(CC)ccc1"
        "smiles": smi
    }
    
    cookie = {
        "csrftoken": "MBAxx5BtO1V9uPWWHMyYimgeP4Jp3kmYdd4nLet3JX5MhDA34DfEwWrYNbcCR4FP"
    }
    # print(data_json['smiles'])
    flag = True
    while flag:
        try:
            re = requests.post(url, data=data_json, cookies=cookie)
            flag = False
        except:
            flag = True
    soup = BeautifulSoup(re.text, "html.parser")
    
    matrix = soup.p.contents[0].split(":")[1]
    matrix = matrix[1:]
    
    return matrix
    # print(soup.p.contents[0])
    
# get_graph_matrix()

# path = "../rdkit_molecule/data/uspto_50.pickle"

# df_data = {
#     "reactants_mol": [], 
#     "products_mol": [],
#     "reaction_type": [],
#     "set": [],
#     "prod_cimg_matrix": [],
#     "reac_cimg_matrix": [],
# }

# df = pd.read_pickle(path)
# # print(df)
# for i in range(df.shape[0]):
#     p_smi = Chem.MolToSmiles(df["products_mol"][i])
#     r_smi = Chem.MolToSmiles(df["reactants_mol"][i])
    
#     df_data["reactants_mol"].append(r_smi)
#     df_data["products_mol"].append(p_smi)
#     df_data["set"].append(df["set"][i])
#     df_data["reaction_type"].append(df["reaction_type"][i])
#     df_data["prod_cimg_matrix"].append(get_graph_matrix(p_smi))
#     df_data["reac_cimg_matrix"].append(get_graph_matrix(r_smi))
#     print(i)
    
# df_new = pd.DataFrame(df_data)
# df_new.to_pickle("../rdkit_molecule/data/uspto_50_graph.pickle")

with open("data/mol_id_6_28.txt", "r") as f:
    lines = f.readlines()
    
    ids = set()
    for line in lines:
        ids.add(line.strip())
        
    print(len(ids))
    
    
with open("data/mol_ids.txt", "w") as f:
    for i in ids:
        f.write(i + "\n")