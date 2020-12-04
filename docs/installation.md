# Installing plycutter

On a recent ubuntu, you can do the following.

        sudo apt-get install libgeos-c1v5 libmpc-dev libspatialindex-dev
        python -m pip install --upgrade pip
        python -m pip install flake8 pytest
        if [ -f requirements.txt ]; then pip install -r requirements.txt; fi

This is copied from the github test action so it should work
ok.

On MacOS, the author uses conda to install the dependencies
but there is no explicit list at the moment.
