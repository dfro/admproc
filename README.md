# AdmProc

AdmProc is Python module for processing measurement files from Admittance program.

## How to use it
```
import admproc

data, freq, area, eps = admproc.read(file_name)
cap, _, voltage, _ = admproc.extract(data, freq, tsel=300, fsel=500000)

```