import os

for module_name in os.listdir("../uptrop"):
    if module_name.endswith(".py"):
        doc_name = module_name.split('.')[0]
        with open("source/" + doc_name + ".rst", 'w') as this_doc:
            this_doc.write("="*len(doc_name) + "\n")
            this_doc.write(doc_name + "\n")
            this_doc.write("="*len(doc_name) + "\n")
            this_doc.write(".. automodule:: uptrop."+doc_name)
            this_doc.write("\n   :members:")


