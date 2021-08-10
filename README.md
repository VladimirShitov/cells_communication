# communiquer
Interface for the several tools that predict cells communication from Single-Cell RNA sequencing data

# Installation for the development
1. Clone the repository
```console
$ git clone https://github.com/VladimirShitov/communiquer.git
```

2. Install [poetry](https://python-poetry.org/). Recommended way:
```console
$ curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python -
```
Alternative way (not recommended):
```console
$ pip install --user poetry
```

3. Install project's python dependencies:
```console
$ poetry install
```

4. Install tools for inferring cells communications. Follow the instructions from the repositories:
* [Cellchat](https://github.com/sqjin/CellChat)
* [CellPhoneDB](https://github.com/Teichlab/cellphonedb)
* [CellCall](https://github.com/ShellyCoder/cellcall)
