## CNRBind

A prediction method by fusing RNA sequence and structure information to identify small molecule-RNA binding sites. 

The RNA id, sequence and labels can be found in ./data_cache/, the codes for RLBind are available in ./predict, the predicted ASA can be found in ./ASA, the PDB files of proteins are saved in ./pdb. Furthermore, the demo and corresponding documentation files can be found in ./demo. See our paper for more details.

Testing each RNA takes approximately 10 seconds, depending on the sequence length.

### Test the model on TE18 (~3min)

```bash
cd ./predict/
python predict_TE18.py
```
### Test the model on RB9 (~3min)
```bash
cd ./predict/
python predict_RB9.py
```
### contact
Kang Xiao: xiaokangneuq@163.com
