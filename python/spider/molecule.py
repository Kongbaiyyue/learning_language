import requests
import json
from bs4 import BeautifulSoup
from rdkit import Chem

import pandas as pd

def get_smi_json():
    # url = "http://www.dcaiku.com/v1/perseus/mol_search"
    # url = "http://www.dcaiku.com/v1/account"
    url = "http://www.dcaiku.com/v1/public/token_refresh"

    headers = {
        "Accept": "application/json, text/plain, */*",
        "Accept-Encoding": "gzip, deflate",
        "Accept-Language": "zh-CN,zh;q=0.9",
        "Authorization":"eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJzdWIiOiJhY2Nlc3NUb2tlbiIsImNyZWF0ZWRBdCI6IjE2ODc5MTQ3NjA0NDEiLCJleHBpcmVkQXQiOiIxNjg3OTIxOTYwNDQxIiwicm9sZSI6IlJPTEVfVVNFUiIsImp0aSI6InBlcnNldXMiLCJlbWFpbCI6Ijk3MjMxMzAwNkBxcS5jb20ifQ.bbCH_izn1j_iDjlNzFzuIcZ-xkMdtR3d3oZ80FMHd2U",
        "Content-Length": "36",
        "Content-Type": "application/json;charset=UTF-8",
        "Host": "www.dcaiku.com",
        "Origin": "http://www.dcaiku.com",
        "Proxy-Connection": "keep-alive",
        "Referer": "http://www.dcaiku.com/",
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/114.0.0.0 Safari/537.36"
    }

    # data_json = [
    #     {"label":"c","type":"title"}
    # ]
    
    data_json = {
        "refreshToken": "eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJzdWIiOiJhY2Nlc3NUb2tlbiIsImNyZWF0ZWRBdCI6IjE2ODc5MTQ3NjA0NDEiLCJleHBpcmVkQXQiOiIxNjg3OTIxOTYwNDQxIiwicm9sZSI6IlJPTEVfVVNFUiIsImp0aSI6InBlcnNldXMiLCJlbWFpbCI6Ijk3MjMxMzAwNkBxcS5jb20ifQ.bbCH_izn1j_iDjlNzFzuIcZ-xkMdtR3d3oZ80FMHd2U"
    }
    
    with open("data/mol_id.txt", "w") as f:
        while True:
            try:
                r = requests.post(url, data=json.dumps(data_json), headers=headers)
                # r = requests.post(url, headers=headers)
                data = json.loads(r.text)
                
                if "error" in data.keys():
                    r = requests.post(url, data=json.dumps(data_json), headers=headers)
                    print(r.text)
                for j in range(len(data["data"]["moleculeItemVos"])):
                    f.write(data["data"]["moleculeItemVos"][j]["id"] + "\n")
                
                data_json["label"] = data["data"]["next"]
                data_json["type"] = "next"
            except:
                continue
    
get_smi_json()

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