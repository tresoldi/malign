coverage run --source src/malign tests/test_malign.py
coverage run --source src/malign -a tests/test_matrix.py
coverage run --source src/malign -a tests/test_demo.py
coverage html
firefox htmlcov/index.html