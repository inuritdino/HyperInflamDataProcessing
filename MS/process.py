from scvi.dataset.dataset10X import Dataset10X
from scvi.dataset.dataset import GeneExpressionDataset

csf_patient_id = ['GSM4104122_MS19270','GSM4104123_MS58637',
                  'GSM4104124_MS71658','GSM4104125_MS49131',
                  'GSM4104126_MS60249','GSM4104127_MS74594']

csf_ctrl_id = ['GSM4104128_PST83775','GSM4104129_PTC32190',
               'GSM4104130_PST95809','GSM4104131_PTC41540',
               'GSM4104132_PST45044','GSM4104133_PTC85037']

dataset = []
for f in csf_patient_id:
    x = Dataset10X(save_path=f)
    dataset.append(x)

CSFMS = GeneExpressionDataset.concat_datasets(*dataset)

dataset = []
for f in csf_ctrl_id:
    x = Dataset10X(save_path=f)
    dataset.append(x)

CSFCTRL = GeneExpressionDataset.concat_datasets(*dataset)
