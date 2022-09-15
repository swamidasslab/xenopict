upload:
	python setup.py build sdist
	twine upload dist/*
	git push --tags
