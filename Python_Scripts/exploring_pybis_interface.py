from aiida_openbis.utils import bisutils
import yaml
import utils

# Connect to openBIS
session = bisutils.log_in(bisurl="openbis", bisuser="admin", bispasswd="changeit")

# Download datasets from openBIS
# ds = session.get_datasets()
# for dataset in ds:
#     dataset.data["dataStore"]["downloadUrl"] = 'https://openbis/'
#     dataset.download(destination = 'work/aiida-openbis/Python_Scripts/openBIS_File')


# Create new datasets in openBIS
s = session.get_sample('20231106213635396-99')
ds_new = session.new_dataset(
    type = 'RAW_DATA',
    sample = '20231106213635396-99',
    files = ['/home/jovyan/work/aiida-openbis/Python_Scripts/openBIS_File/20231106155424329-61/original/SI Edge-extended zigzag graphene nanoribbons.pdf']
)
ds_new.save()