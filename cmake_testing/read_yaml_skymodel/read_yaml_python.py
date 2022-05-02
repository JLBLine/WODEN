import yaml

with open(r'point_powerlaw.yaml') as file:
    documents = yaml.full_load(file)

    for item, doc in documents.items():
        # print(item, ":", doc)

        print(type(doc[0]))
