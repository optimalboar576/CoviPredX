import os
import pandas as pd
import subprocess
import tkinter as tk
from tkinter import filedialog, messagebox
from PIL import Image, ImageTk  # For adding images
from rdkit import Chem
from rdkit.Chem import AllChem
import xgboost as xgb

class BioactivityScreeningApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Bioactivity Screening Tool")
        self.root.geometry("600x500")

        # Use relative path for the logo image directly
        logo_path = "data/logo1.jpg"
        self.logo_image = Image.open(logo_path)
        self.logo_image = self.logo_image.resize((450, 300), Image.LANCZOS)
        self.logo = ImageTk.PhotoImage(self.logo_image)
        self.logo_label = tk.Label(root, image=self.logo)
        self.logo_label.pack(pady=10)

        # Create a label and button to upload a file
        self.upload_button = tk.Button(root, text="Upload CSV File", font=("Helvetica", 12), command=self.upload_file)
        self.upload_button.pack(pady=10)

        # Label to show file upload status
        self.file_label = tk.Label(root, text="", font=("Helvetica", 10))
        self.file_label.pack(pady=5)

        # Button to calculate descriptors and make predictions
        self.predict_button = tk.Button(root, text="Calculate Descriptors & Predict", font=("Helvetica", 12), command=self.process_file)
        self.predict_button.pack(pady=10)

        # Label to show processing status
        self.status_label = tk.Label(root, text="", font=("Helvetica", 10))
        self.status_label.pack(pady=5)

        self.filepath = None
        self.predictions_df = None

    def upload_file(self):
        self.filepath = filedialog.askopenfilename(title="Select CSV file", filetypes=(("CSV Files", "*.csv"),))
        if self.filepath:
            self.file_label.config(text=f"File selected: {os.path.basename(self.filepath)}")
        else:
            self.file_label.config(text="No file selected")

    def calculate_morgan_fingerprints(self, smiles_data):
        morgan_fps = []
        for smi in smiles_data["SMILES"]:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
                morgan_fps.append(list(fp))
            else:
                morgan_fps.append([None] * 1024)

        morgan_df = pd.DataFrame(morgan_fps, columns=[f'fp_{i}' for i in range(1024)])
        morgan_df["ID"] = smiles_data["ID"].values
        return morgan_df

    def desc_calc(self, smiles_file):
        smiles_data = pd.read_csv(smiles_file, sep=',', header=None, names=["SMILES", "ID"])
        smiles_data.to_csv('molecule.smi', sep='\t', header=False, index=False)

        # Use relative path for PaDEL descriptor jar and descriptor file
        padel_jar_path = "data/PaDEL-Descriptor/PaDEL-Descriptor.jar"
        descriptor_type = "data/PaDEL-Descriptor/PubchemFingerprinter.xml"

        bashCommand = f"java -Xms2G -Xmx2G -Djava.awt.headless=true -jar {padel_jar_path} -removesalt -standardizenitro -fingerprints -descriptortypes {descriptor_type} -dir . -file padel_descriptors_output.csv"
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        process.communicate()

        # Read PaDEL descriptors
        padel_descriptors = pd.read_csv('padel_descriptors_output.csv')
        morgan_fps = self.calculate_morgan_fingerprints(smiles_data)

        padel_descriptors['Name'] = padel_descriptors['Name'].astype(str)
        morgan_fps['ID'] = morgan_fps['ID'].astype(str)

        combined_df = pd.merge(padel_descriptors, morgan_fps, left_on="Name", right_on="ID")
        top_n_features_df = pd.read_csv('data/top_n_features.csv')
        top_n_feature_names_list = top_n_features_df['Feature'].tolist()

        important_columns = ['Name'] + top_n_feature_names_list
        filtered_combined_df = combined_df[important_columns]

        filtered_combined_df.to_csv('filtered_combined_descriptors.csv', index=False)
        return filtered_combined_df

    def prediction(self, input_data):
        molecule_name = input_data["Name"]
        input_data = input_data.drop(columns=["Name"])

        # Use relative path for the model
        model_path = "data/Replicase_xgb_model.json"
        load_model = xgb.Booster()
        load_model.load_model(model_path)

        dmatrix_input = xgb.DMatrix(input_data)
        prediction = load_model.predict(dmatrix_input)

        prediction_output = pd.Series(prediction, name='pIC50')
        df = pd.concat([molecule_name.reset_index(drop=True), prediction_output], axis=1)
        return df

    def process_file(self):
        if not self.filepath:
            messagebox.showerror("Error", "Please upload a file first")
            return

        self.status_label.config(text="Processing...")

        try:
            # Process the file to calculate descriptors and make predictions
            desc_data = self.desc_calc(self.filepath)
            self.predictions_df = self.prediction(desc_data)

            # Save the prediction file in the same directory as the input file
            output_file = os.path.join(os.path.dirname(self.filepath), 'predictions.csv')
            self.predictions_df.to_csv(output_file, index=False)

            messagebox.showinfo("Success", f"Prediction complete! Results saved to {output_file}")
            self.status_label.config(text="Prediction complete!")
        except Exception as e:
            messagebox.showerror("Error", str(e))
            self.status_label.config(text="Error in processing")
        finally:
            # Remove temporary files created during processing
            try:
                os.remove("filtered_combined_descriptors.csv")
                os.remove("molecule.smi")
                os.remove("padel_descriptors_output.csv")
                print("Temporary files removed successfully.")
            except Exception as e:
                print(f"Error removing temporary files: {e}")

# Create the root window
root = tk.Tk()
app = BioactivityScreeningApp(root)

# Run the application
root.mainloop()
