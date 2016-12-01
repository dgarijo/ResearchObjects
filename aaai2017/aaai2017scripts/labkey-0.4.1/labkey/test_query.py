from labkey.utils import create_server_context
from labkey.query import select_rows, QueryFilter

print("Create a server context")
server_context = create_server_context('idri.labkey.com', 'Formulations')

ssr = select_rows(server_context, 'Samples', 'Formulations')
print("select_rows: There are " + str(len(ssr['rows'])) + " rows.")

ssr = select_rows(server_context, 'Samples', 'Formulations', filter_array=[
    QueryFilter('Batch', 'QF', filter_type=QueryFilter.Types.CONTAINS)
])
print("select_rows: There are " + str(len(ssr['rows'])) + " filtered rows.")