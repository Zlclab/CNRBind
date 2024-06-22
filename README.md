## CNRBind: Small molecule-RNA binding sites recognition via site significance from nucleotide and complex network information. 

The RNA id, sequence and labels can be found in ./data_cache.
The codes for CNRBind are available in ./predict.
The predicted ASA can be found in ./ASA.
The PDB files are saved in ./pdb.
The predicted PDB files are saved in ./predicted_pdb.

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
