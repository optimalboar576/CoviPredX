# CoviPredX - Bioactivity Screening Tool

Welcome to **CoviPredX**, a bioactivity screening tool designed to predict the pIC50 values for COVID-19-related molecules using machine learning models and molecular fingerprints. This tool provides a user-friendly interface for researchers and scientists to predict bioactivity, offering accurate and fast predictions. A web application of CoviPredX is availabe on https://covipredx.bicpu.edu.in

## Why Use CoviPredX?

- **Accurate Predictions**: Built with cutting-edge machine learning algorithms like XGBoost, CoviPredX provides reliable pIC50 predictions.
- **Easy to Use**: With a graphical user interface (GUI) powered by Tkinter, the tool allows you to upload CSV files and get predictions with just a few clicks.
- **Versatile Input**: Accepts CSV files containing SMILES strings, making it suitable for a wide range of molecular input data.
- **Comprehensive Descriptor Calculation**: Automatically calculates Morgan fingerprints and PaDEL descriptors.
- **No Coding Required**: Aimed at researchers without programming expertise, CoviPredX allows bioactivity predictions through a simple interface.

## Download (For Windows/Ubuntu)

To get started, download the latest version of **CoviPredX** from the [Releases](https://github.com/your-repo/releases) section on GitHub. It is also available in the Downloads section of the web-server.

## OR

## Clone the Repository

To clone the repository and navigate to the project directory:

```bash
git clone https://github.com/your-repo/CoviPredX.git
cd CoviPredX
```

## Installation Instructions

### For Ubuntu

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
