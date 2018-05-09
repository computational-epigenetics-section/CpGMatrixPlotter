# CpGMatrixPlotter
package to plot CpG matrix data as circles spaced proportionally apart

## Plot example
![Plot example](cpgPlotter/examples/pcdh_region_chr18_37449400.png)

## Requirements
* Python3+
* Numpy
* Matplotlib

## Install
The package is not yet uploaded to PyPi so you must clone/download the repo first
1. Git clone or download and unzip the repo
2. Activate your virtualenv (Optional)
3. Navigate to the directory containing the `setup.py` file
4. Execute `pip install .` (note the period at the end)

## Usage
```python
from cpgPlotter import CpgMatrixPlotter
data = np.array([[1,1,1,0],
                 [1,1,0,0]])
locations = np.array([100, 110, 155, 190])

plotter = CpgMatrixPlotter()
plotter.plotCpgMatrix(data, locations)

```
