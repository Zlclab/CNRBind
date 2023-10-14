## CNRBind

A prediction method by fusing RNA sequence and structure information to identify small molecule-RNA binding sites. 

The RNA id, sequence and labels can be found in ./data_cache/, the codes for RLBind are available in ./predict. In addition, the demo and corresponding documentation files can be found in ./demo.

### Test the model on TE18 (~3min)

```bash
cd ./predict/
python predict_TE18.py
```
### Re-training your own model for the new dataset
```bash
cd ./predict/
python predict_RB9.py
```
### contact
Kang Xiao: xiaokangneuq@163.com
