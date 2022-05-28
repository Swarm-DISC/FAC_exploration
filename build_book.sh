# Compiles the Jupyter Book from .ipynb & .py
#  1. Re-introduce a standardised kernelspec to nb metadata (cleared by pre-commit)
#  2. Convert .py files to .md, adding annotations in each, and referenced in _toc.yml
#  3. Use jupyter-book to execute notebooks and create the html

# Fix the notebooks' metadata
python adjust_nb_metadata.py

# Create the modified .py files as .md to include in book
echo "Creating temporary .md files in _prebuild directory"
mkdir -p _pre_build
rm -rf _pre_build/*
annotation="**NB: Do not use the download button above. Go to the repository on GitHub to access the files.**"
for f in "funcs_fac.py" "plot_3sat_conjunction.py" "plot_and_save_dual_sat.py" "plot_and_save_single_sat_MVA.py" "plot_and_save_single_sat.py" "plot_and_save_three_sat.py" "plot_qi.py" "save_qi.py"
do
printf '%s\n\n%s\n\n```python\n%s\n```' "# $f" "$annotation" "$(cat notebooks/$f)" > _pre_build/$f.md
echo "Converted $f to $f.md"
done

# Build the book
jupyter-book build .
