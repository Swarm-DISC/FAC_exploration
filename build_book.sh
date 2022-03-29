# Compiles the Jupyter Book from .ipynb & .py
#  1. Clear notebook outputs, including kernel metadata
#     (removing the kernelspec helps ensure that the current kernel is used instead)
#  2. Convert .py files to .md, referenced in _toc.yml
#  3. Build the book, caching the notebooks' executed states
pre-commit run --all
# jupyter-book clean . --all
echo "Creating temporary .md files in _prebuild directory"
mkdir -p _pre_build
rm -rf _pre_build/*
annotation="**NB: Do not use the download button above. Go to the repository on GitHub to access the files.**"
for f in "funcs_fac.py" "plot_3sat_conjunction.py" "plot_and_save_dual_sat.py" "plot_and_save_single_sat_MVA.py" "plot_and_save_single_sat.py" "plot_and_save_three_sat.py" "plot_qi.py" "save_qi.py"
do
printf '%s\n\n%s\n\n```python\n%s\n```' "# $f" "$annotation" "$(cat notebooks/$f)" > _pre_build/$f.md
echo "Converted $f to $f.md"
done
jupyter-book build .
