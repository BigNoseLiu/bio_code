#/root/miniconda3/bin/python3 manage.py  startapp cancer_database
#/root/miniconda3/bin/python3 manage.py  startapp inher_data
/root/miniconda3/bin/python3 manage.py makemigrations
/root/miniconda3/bin/python3 manage.py migrate

/root/miniconda3/bin/python3 manage.py runserver 10.10.9.22:6001
