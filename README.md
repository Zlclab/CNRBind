## CNRBind

A prediction method by fusing RNA sequence and structure information to identify small molecule-RNA binding sites. 

The RNA id, sequence and labels can be found in ./data_cache
The codes for RLBind are available in ./predict
The predicted ASA can be found in ./ASA
The PDB files of proteins are saved in ./pdb

Testing each RNA takes approximately 10 seconds, depending on the sequence length.

### Test the model on TE18 (~3min)

```bash
cd ./predict/
python predict.py TE18
```
### Test the model on RB9 (~3min)
```bash
cd ./predict/
python predict.py RB9
```
### contact
Kang Xiao: xiaokangneuq@163.com
