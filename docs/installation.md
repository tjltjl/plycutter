# Installing plycutter

On a recent ubuntu, you can do the following.

        sudo apt-get install libgeos-c1v5 libmpc-dev libspatialindex-dev
        python -m pip install --upgrade pip
        python -m pip install flake8 pytest
        if [ -f requirements.txt ]; then pip install -r requirements.txt; fi

This is copied from the github test action so it should work
ok.

On MacOS (and probably other OSes), you can use conda to install the dependencies using

        conda env update --file environment.yml
        
in the main plycutter directory. If you like, you can create a new conda environment
for plycutter before doing this.
