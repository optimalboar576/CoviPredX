# CoviPredX - Bioactivity Screening Tool

Welcome to **CoviPredX**, a bioactivity screening tool designed to predict the pIC50 values for COVID-19-related molecules using machine learning models and molecular fingerprints. This tool provides a user-friendly interface for researchers and scientists to predict bioactivity, offering accurate and fast predictions. A web application of CoviPredX is availabe on http://covipredx.bicpu.edu.in

## Why Use CoviPredX?

- **Accurate Predictions**: Built with cutting-edge machine learning algorithms like XGBoost, CoviPredX provides reliable pIC50 predictions.
- **Easy to Use**: With a graphical user interface (GUI) powered by Tkinter, the tool allows you to upload CSV files and get predictions with just a few clicks.
- **Versatile Input**: Accepts CSV files containing SMILES strings, making it suitable for a wide range of molecular input data.
- **Comprehensive Descriptor Calculation**: Automatically calculates Morgan fingerprints and PaDEL descriptors.
- **No Coding Required**: Aimed at researchers without programming expertise, CoviPredX allows bioactivity predictions through a simple interface.

## Download (For Windows/Ubuntu)

To get started, download the latest version of **CoviPredX** from the [Releases](http://github.com/optimalboar576/CoviPredX/releases) section on GitHub. It is also available in the Downloads section of the [web-server](http://covipredx.bicpu.edu.in).

## OR

## Clone the Repository

To clone the repository and navigate to the project directory:

```bash
git clone https://github.com/optimalboar576/CoviPredX.git
cd CoviPredX
```

## Installation Instructions

## 1. For Ubuntu

### First unzip and navigate to the CoviPredX directory
```bash
unzip CoviPredX.zip
cd CoviPredX
```
### Create an environment using the provided environment.yml file....
```bash
conda env create -f environment.yml
```
### Activate the environment
```bash
conda activate covipredx
```
### Install the tool
```bash
pip install covipredx-1.0.0-py3-none-any.whl
```
### Check if installation worked
```bash
covipredx --help
```
### USAGE
```bash
usage: covipredx [-h] input_file output_file
```

## 2. For Windows
It is a portable windows application, which needs no installation. Simply unzip the CoviPredX.gz and click on the CoviPredX.bat batch file to run the application. Make sure you have installed the following prerequisites before running the application--

- **Python 3.10 or higher**: Download from [here](https://www.python.org/downloads/).
- **Java Runtime Environment (JRE)**: Required for PaDEL descriptor calculations. Download from [here](https://www.oracle.com/java/technologies/javase-downloads.html).
- **pillow**
- **pandas**
- **numpy**
- **xgboost**
- **rdkit**
