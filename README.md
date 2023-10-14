# CNRBind
A prediction method by fusing RNA sequence and structure information to identify small molecule-RNA binding sites. 
The RNA id, sequence and labels are available in ./data_cashe. 
The PDB files of RNAs are availavle in ./pdb. 
The codes for CNRBind are available in ./src. 
The demo and corresponding documentation files can be found in ./demo. 
See our paper for more details.

Testing each RNA takes approximately 10 seconds, depending on the sequence length.

## Test the model on TE18 (~3min)
```hello world``` 
'<cd ./predict/>'
python predict_TE18.py  

## Test the model on RB9 (~3min)
cd ./predict/
python predict_RB9.py  


Kang Xiao: xiaokangneuq@163.com
