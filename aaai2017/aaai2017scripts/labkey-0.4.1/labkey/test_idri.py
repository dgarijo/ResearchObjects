# Sample to get a batch

from utils import create_server_context
from experiment import load_batch

print("Create a server context")
project_name = 'Formulations'  # Project folder name
server_context = create_server_context('idri.labkey.com', project_name, use_ssl=True)

print("Load an Assay batch from the server")
assay_id = 97  # provide one from your server
batch_id = 4777  # provide one from your server
run_group = load_batch(assay_id, batch_id, server_context)

if run_group is not None:
    print("Batch Id: " + str(run_group.id))
    print("Created By: " + run_group.created_by)

