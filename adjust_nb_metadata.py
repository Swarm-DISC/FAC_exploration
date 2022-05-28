"""Rewrite the metadata of notebooks

Standardises the kernelspec so that
1. jupyter-book uses the right kernel in CI for execution
2. jupyter-book identifies them as notebooks,
    so that the launcher buttons are added
"""

from glob import glob
import nbformat

# Load each notebook
notebooks = glob("./notebooks/*.ipynb", recursive=True)
for ipath in notebooks:
    print("Rewriting metadata for", ipath)
    ntbk = nbformat.read(ipath, nbformat.NO_CONVERT)
    # Insert normalised kernelspec
    ntbk["metadata"]["kernelspec"] = {
        "display_name": "Python 3",
        "language": "python",
        "name": "python3"
    }
    # Overwrite notebook with the corrected content
    nbformat.write(ntbk, ipath)
