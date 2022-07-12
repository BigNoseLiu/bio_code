#sudo pip3 install Django
sudo /root/miniconda3/bin/pip install Django
/root/miniconda3/bin/pip install django-import-export
/root/miniconda3/bin/django-admin startproject bio_web

ln -s site-packages/django/contrib/admin/static/admin/css/base.css /root/miniconda3/lib/python3.7/site-packages/django/contrib/admin/static/admin/css/base.css
ln -s site-packages/django/contrib/admin/static/admin/css/forms.css /root/miniconda3/lib/python3.7/site-packages/django/contrib/admin/static/admin/css/forms.css
