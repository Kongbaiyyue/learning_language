import requests
import json


ccf_A_list = [
    "https://icml.cc/static/virtual/data/icml-2021-orals-posters.json",
    "https://nips.cc/static/virtual/data/neurips-2021-orals-posters.json",
    "https://iclr.cc/static/virtual/data/iclr-2021-orals-posters.json",

]

AAAI = "https://aaai-2022.virtualchair.net/papers.json"
IJCAI = "https://www.ijcai.org/proceedings/2022/"


def get_paper(url, name):
    content = requests.get(url).content
    result = json.loads(content)
    path = f"./{name}.json"
    with open(path,"w") as f :
        json.dump(result,f)
    # print(content)
    return path


url = "https://icml.cc/static/virtual/data/icml-2021-orals-posters.json"
name = "icml2021"
path = get_paper(url, name)

def filter_paper_name(path, name):
    with open(path, "r") as f1 , open(f"./{name}_paper_titles.txt", "w", encoding="utf-8") as f2:
        all_data = json.load(f1)

        paper_titles = set()
        for paper in all_data["results"]:
            paper_title = paper["name"]
            paper_titles.add(paper_title)
        # print(paper_titles)

        for title in paper_titles:
            f2.write(title + "\n")


filter_paper_name(path, name)