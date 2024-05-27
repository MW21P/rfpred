![Project Logo](assets/banner.png)

![Coverage Status](assets/coverage-badge.svg)

<h1 align="center">
<i><b>rfpred 🧪 </b></i>
</h1>

<br>


🔍 Package that predicts the Rf values for silica thin-layer-chromatographies (TLCs) in synthesis labs based on the chemical compound and the solvent system.

## 🔥 Usage

After installing, you should create a new file and activate the environment where you installed ***rfpred***. This then allows you to use our prediction mode on a local host by copying this code into your new file:

```python
import rfpred

# One line to rule them all
rfpred.App.run()
```

## 👩‍💻 Installation

Create a new environment, you may also give the environment a different name and activate this new environment.

```
conda create -n rfpred python=3.10
```
```
conda activate rfpred
```
Then simply pip install the package by copying this into your command line.
```
(rfpred) $ pip install "git+https://github.com/MW21P/rfpred.git"
```
Then create a new file and proceed how described in the Usage section.

### Run tests and coverage

```
(conda_env) $ pip install tox
(conda_env) $ tox
```

### 📖 Authors
Milena Wiegand: https://github.com/MW21P

Matthias Galka: https://github.com/MGalka66

This project was carried out as part of EPFL's ***super cool*** "Practical programming in Chemistry" course.

