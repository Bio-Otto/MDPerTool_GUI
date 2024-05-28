::pip install . -vv
pip install pillow==9.0.0
python setup.py install --single-version-externally-managed --record=record.txt

::"%PYTHON%" setup.py install
::if errorlevel 1 exit 1