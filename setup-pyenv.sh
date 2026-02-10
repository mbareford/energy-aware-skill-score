PRFX=${HOME}/utils/pyenv
cd ${PRFX}


https://repo.anaconda.com/miniconda/Miniconda3-py313_25.11.1-1-Linux-x86_64.sh

PYTHON_LABEL=py313
MINICONDA_TAG=miniconda
MINICONDA_LABEL=${MINICONDA_TAG}3
MINICONDA_VERSION=25.11.1-1
MINICONDA_ROOT=${PRFX}/${PYTHON_LABEL}

cd ${PRFX}

MINICONDA_TITLE=${MINICONDA_LABEL^}
MINICONDA_BASH_SCRIPT=${MINICONDA_TITLE}-${PYTHON_LABEL}_${MINICONDA_VERSION}-Linux-x86_64.sh

wget https://repo.anaconda.com/${MINICONDA_TAG}/${MINICONDA_BASH_SCRIPT}
chmod 700 ${MINICONDA_BASH_SCRIPT}
unset PYTHONPATH
bash ${MINICONDA_BASH_SCRIPT} -b -f -p ${MINICONDA_ROOT}
rm ${MINICONDA_BASH_SCRIPT}
cd ${MINICONDA_ROOT}

PATH=${MINICONDA_ROOT}/bin:${PATH}
conda init --dry-run --verbose > activate.sh
conda_env_start=`grep -n "# >>> conda initialize >>>" activate.sh | cut -d':' -f 1`
conda_env_stop=`grep -n "# <<< conda initialize <<<" activate.sh | cut -d':' -f 1`

echo "sed -n '${conda_env_start},${conda_env_stop}p' activate.sh > activate2.sh" > sed.sh
echo "sed 's/^.//' activate2.sh > activate.sh" >> sed.sh
echo "rm activate2.sh" >> sed.sh
. ./sed.sh
rm ./sed.sh

. ${MINICONDA_ROOT}/activate.sh

conda update -y -n root --all

export PS1="(pyenv) [\u@\h \W]\$ "

pip install --upgrade pip
pip install numpy
pip install scipy
pip install pandas
pip install matplotlib
pip install netCDF4
pip install xarray


mkdir ./repos
cd ./repos
git clone https://github.com/w-k-jones/get_jasmin_era5.git
cd get_jasmin_era5
python setup.py install


conda deactivate
export PS1="[\u@\h \W]\$ "
