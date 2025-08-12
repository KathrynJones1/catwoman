rm -rf build dist *.egg-info
python -m build
twine check dist/*
twine upload dist/*
