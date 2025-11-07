

import os
import sys
import subprocess
import pkg_resources
import cmath  # explicitly import to ensure bundlers include this stdlib module
import tkinter as tk
from tkinter import filedialog, messagebox

def install_dependencies(missing_packages):
    """Install missing packages using pip."""
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install"] + missing_packages)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error installing packages: {e}")
        return False

def check_and_install_dependencies():
    """Check if required packages are installed and offer to install if missing."""
    required_packages = {
        'pandas': 'pandas>=2.0.0',
        'numpy': 'numpy>=1.24.0',
        'matplotlib': 'matplotlib>=3.7.0',
        'openpyxl': 'openpyxl>=3.1.0'  # Required for reading xlsx files
    }
    
    missing = []
    for package, spec in required_packages.items():
        try:
            pkg_resources.require(spec)
        except (pkg_resources.VersionConflict, pkg_resources.DistributionNotFound):
            missing.append(spec)
    
    if missing:
        print("Missing required packages:", missing)
        root = tk.Tk()
        root.withdraw()  # Hide the main window
        
        if messagebox.askyesno("Missing Dependencies", 
            "Some required packages are missing. Would you like to install them now?\n\n" +
            "Missing packages:\n" + "\n".join(missing)):
            
            if install_dependencies(missing):
                messagebox.showinfo("Success", "Dependencies installed successfully. Please restart the application.")
            else:
                messagebox.showerror("Error", "Failed to install dependencies. Please install them manually:\n\npip install " + " ".join(missing))
            
            root.destroy()
            sys.exit(1)
        else:
            messagebox.showwarning("Warning", "The application may not work correctly without required packages.")
            root.destroy()
            sys.exit(1)

# Check and install dependencies before any other imports
if __name__ == '__main__':
    # If running as a frozen app (py2app/PyInstaller), dependencies should be bundled.
    if getattr(sys, 'frozen', False):
        try:
            import matplotlib.pyplot as plt
            import pandas as pd
            import numpy as np
        except Exception as e:
            # If a package is missing from the bundle, inform the user and instruct to rebuild.
            tmp_root = tk.Tk()
            tmp_root.withdraw()
            messagebox.showerror(
                "Missing bundled dependency",
                f"A required package is missing from the application bundle: {e}\n\n"
                "Rebuild the app and ensure these packages are included (see setup.py OPTIONS['packages'])."
            )
            sys.exit(1)
    else:
        # In development mode, offer to install missing packages via pip
        check_and_install_dependencies()
        try:
            import matplotlib.pyplot as plt
            import pandas as pd
            import numpy as np
        except ImportError as e:
            tmp_root = tk.Tk()
            tmp_root.withdraw()
            messagebox.showerror("Import Error", f"Failed to import required packages: {e}\n\nPlease install dependencies and restart.")
            sys.exit(1)

    # Start the main application
    root = tk.Tk()
    root.geometry("800x400")
    root.title("ATP Production Calculator")
    label = tk.Label(root, text="ATP Production Calculator")
    label.grid(column=0, row=0)


# add your GUI components and logic here

def atp_prod_calc(file_path, res_dir,oligo_inj, rot_inj, run_acid=1):
    print("Processing file:", file_path)
    label.configure(text="Processing file: " + label_file_explorer.cget("text").replace("File Opened: ",""))
        
    # Read assay configuration, rate, and raw data sheets from the excel file
    assay_conf = pd.read_excel(file_path, sheet_name="Assay Configuration")
    rate = pd.read_excel(file_path, sheet_name="Rate")
    raw = pd.read_excel(file_path, sheet_name="Raw")

    # Extract plate layout and groups from assay configuration
    plate_idx = assay_conf.index[assay_conf.iloc[:,0] == "Group Layout"].to_list()
    plate_layout = assay_conf.iloc[plate_idx[0]+1:plate_idx[0]+9,2:14]
    plate_layout.columns = np.arange(1, 13).tolist()

    groups=plate_layout.melt().value.unique().tolist()
    print(len(groups), "Experimental groups found")
    
    # extract protocol information and oligomycin and rotenone/antimycin A injection times
    prot_idx = assay_conf.index[assay_conf.iloc[:,0] == "Protocol Summary"].to_list()
    protocol = assay_conf.iloc[prot_idx[0]+1:prot_idx[0]+8, 2:7]

    cycle_counts = protocol.iloc[5, :].astype(str).str[-1].astype(int)
    injections = protocol.iloc[0,:].tolist()
    
    # Store both the injection names and their indices
    oligomycin = injections[oligo_inj]  # Full injection name
    rot_aa = injections[rot_inj]  # Full injection name
    print("The following will be used as the oligomycin injection:", oligomycin)
    print("The following will be used as the rotenone/antimycin A injection:", rot_aa)
    
    protocol_df = pd.DataFrame(zip(injections, cycle_counts), columns=["Injection", "Cycle Count"]).set_index("Injection")
    protocol_df["cycle_end"] = protocol_df["Cycle Count"].cumsum()
    protocol_df["cycle_start"] = protocol_df["cycle_end"] - protocol_df["Cycle Count"] + 1

    if run_acid:
        print("Running buffering power calculation...")
        # check groups for sulfuric acid and hydrochloric acid controls
        sulfuric_acid_ctrl = [g for g in groups if "sul" in str(g).lower() or "h2so4" in str(g).lower()]
        hydrochloric_acid_ctrl = [g for g in groups if "hcl" in str(g).lower()]

        if sulfuric_acid_ctrl != []:
            print(f"The following groups appear to be sulfuric acid controls: {sulfuric_acid_ctrl}")
        if hydrochloric_acid_ctrl != []:
            print(f"The following groups appear to be hydrochloric acid controls: {hydrochloric_acid_ctrl}")
        if sulfuric_acid_ctrl == [] and hydrochloric_acid_ctrl == []:
            sulfuric_acid_ctrl = [g for g in groups if "buff" in str(g).lower() or "acid" in str(g).lower()]
            print(f"The following groups will be used as sulfuric acid controls: {sulfuric_acid_ctrl}")
        else:
            sulfuric_acid_ctrl = []
            run_acid = 0
            print("No sulfuric acid controls detected. Standard buffering power of 0.28939800915432956 will be used.")
    
        # calculating buffering power. For each buffering power well calculate the average pH reading from the final 5 readings prior to any injection
        raw_df = raw.loc[raw.Group == sulfuric_acid_ctrl[0], ["pH", "Measurement", "Tick", "Well", "Group"]]
        basal_pH_cycle = protocol_df.iloc[0, 1]
        print("collecting basal pH reading from cycle:", basal_pH_cycle)
    
        # Calculate average pH for the last 5 measurements per well at the basal cycle
        basal_pH_by_well = (raw_df.loc[raw_df.Measurement == basal_pH_cycle]
                            .groupby('Well')
                            .apply(lambda x: x.tail(5)['pH'].mean())
                            .to_dict())
        # Calculate the change in pH compared with the basal pH for each buffering power well using the final 5 readings of each final measurement cycle per injection
        # Store results in a dataframe
        results_list = []
        valid_acid = 1
        try:
            for injection, cycle_end in protocol_df['cycle_end'].items():
                # Calculate mean pH for the last 5 measurements of this cycle for each well
                injection_pH_by_well = (raw_df.loc[raw_df.Measurement == cycle_end]
                                    .groupby('Well')
                                    .apply(lambda x: x.tail(5)['pH'].mean())
                                    .to_dict())
                
                # Create DataFrame rows for this injection
                for well, ph in injection_pH_by_well.items():
                    delta_pH = ph - basal_pH_by_well[well]
                    results_list.append({
                        'Injection': injection,
                        'Well': well,
                        'pH': ph,
                        'Basal_pH': basal_pH_by_well[well],
                        'ΔpH': delta_pH
                    })

            # Convert list of dictionaries to DataFrame
            results_df = pd.DataFrame(results_list)

            # Display the results
            results_df.to_csv(os.path.join(res_dir, f"delta_pH_results.csv"), index=False)
            print(f"ΔpH results saved to CSV delta_pH_results.csv in {res_dir}.")
        except KeyError as e:
            valid_acid = 0
            print("Error calculating ΔpH results:", e)
    # Plot the change in pH against the amount of acid added per well.
    # The standard protocol adds 62.5 nmol H+ per injection. We set baseline (injections[0]) to 0,
    # and for subsequent injections multiply 62.5 by the injection index (0-based).

    if valid_acid and run_acid:
        # Build a mapping from injection name to 0-based injection index
        injection_index = {inj: idx for idx, inj in enumerate(injections)}

        # Compute acid added per injection for each row: 62.5 nmol * injection_index (baseline -> 0)
        results_df['Acid_Added_nmol'] = results_df['Injection'].map(injection_index).fillna(0).astype(int) * 62.5

        # convert nmol H+ to nmol/μL: volume begins at 150 μL and increases by 25 μL per injection
        results_df['Volume_uL'] = 150 + results_df['Injection'].map(injection_index).fillna(0).astype(int) * 25
        results_df['Acid_Added_nmol_per_uL'] = 2.28 * (results_df['Acid_Added_nmol'] / results_df['Volume_uL'])

        plt.figure(figsize=(8, 6))
        for well in results_df['Well'].unique():
            well_data = results_df[results_df['Well'] == well]
            plt.plot(well_data['Acid_Added_nmol_per_uL'], well_data['ΔpH'], marker='o', label=f'Well {well}')
        plt.xlabel('Acid Added (nmol H+/μL)')
        plt.ylabel('Change in pH (ΔpH)')
        plt.title('Change in pH vs. Acid Added per Well')
        plt.legend()
        plt.grid()
        plt.savefig(os.path.join(res_dir, "buffering_power_plot.png"))
        # Linear regression to extract slope (buffering power) for each well
        # First, we need to check if the cycle count for the first injection is longer than standard (i.e., >3)
        # If yes, we will ignore the baseline measurement cycles for buffering power calculation
        buffering_power = {}
        for well in results_df['Well'].unique():
            well_data = results_df[results_df['Well'] == well]
            # Exclude baseline cycle if first injection cycle count > 3
            if protocol_df.iloc[1, 0] > 3:
                well_data = well_data[well_data['Injection'] != injections[0]]
            # Perform linear regression (1st degree polynomial fit)
            slope, intercept = np.polyfit(well_data['Acid_Added_nmol_per_uL'], well_data['ΔpH'], 1)
            buffering_power[well] = abs(slope)  # slope is ΔpH per nmol H+. This is equal to ΔmpH per pmol H+.
    # Convert ECAR from mpH/min to pmol H+/min.
    # Experimental ECAR and OCR data are in the rate dataframe
    # Calculate the average of the last two ECAR measurements prior to oligomycin injection for each experimental well ("Group" is not buffering power or background)
    #     Divide this by the average buffering power to yield the proton production rate in pmol H+/min.

    # average buffering power across all buffering power wells
        avg_buffering_power = np.mean(list(buffering_power.values()))
        print("Average buffering power (ΔpH per nmol H+):", avg_buffering_power)
    else:
        avg_buffering_power = 0.28939800915432956  # default value if no valid acid control is present
        print("Using default average buffering power (ΔpH per nmol H+):", avg_buffering_power)

    ecar_pmol = {}
    for well in rate['Well'].unique():
        well_data = rate[rate['Well'] == well]
        # Get the last two ECAR measurements prior to oligomycin injection
        ecar_last_two = well_data[well_data['Measurement'] <= protocol_df.loc[oligomycin, "cycle_start"]-1]['ECAR'][-2:]
        ecar_pmol[well] = ecar_last_two.mean() / avg_buffering_power  # Average of the last two measurements

    print("\nProton production rates (pmol H+/min) by well calculated from ECAR using cycle numbers:", protocol_df.loc[oligomycin, "cycle_start"]-2, protocol_df.loc[oligomycin, "cycle_start"]-1)

    # We need to calculate the untreated oxygen consumption rate (OCR) for each well using the average of the last two OCR measurements prior to oligomycin injection
    ocr_pmol = {}
    for well in rate['Well'].unique():
        well_data = rate[rate['Well'] == well]
        # Get the last two OCR measurements prior to oligomycin injection
        ocr_last_two = well_data[well_data['Measurement'] <= protocol_df.loc[oligomycin, "cycle_start"]-1]['OCR'][-2:]
        ocr_pmol[well] = ocr_last_two.mean()  # Average of the last two measurements

    # Next we need the OCR for each well in the last two measurements after oligomycin injection
    ocr_oligo_pmol = {}
    for well in rate['Well'].unique():
        well_data = rate[rate['Well'] == well]
        # Get the last two OCR measurements after oligomycin injection
        ocr_oligo_last_two = well_data[well_data['Measurement'] <= protocol_df.loc[oligomycin, "cycle_end"]]['OCR'][-2:]
        ocr_oligo_pmol[well] = ocr_oligo_last_two.mean()  # Average of the last two measurements

    # Next we need the OCR for each well in the last two measurements after rotenone/antimycin A injection
    ocr_rot_pmol = {}
    for well in rate['Well'].unique():
        well_data = rate[rate['Well'] == well]
        # Get the last two OCR measurements after rotenone/antimycin A injection
        ocr_rot_last_two = well_data[well_data['Measurement'] <= protocol_df.loc[rot_aa, "cycle_end"]]['OCR'][-2:]
        ocr_rot_pmol[well] = ocr_rot_last_two.mean()  # Average of the last two measurements

    # Now we can calculate basal OCR by subtracting OCR_rot from the untreated OCR (that is, prior to oligomycin)
    basal_ocr = {}
    for well in rate['Well'].unique():
        basal_ocr[well] = ocr_pmol[well] - ocr_rot_pmol[well]

    # Similarly, we can calculated ATP-linked OCR by subtracting OCR_oligo from the untreated OCR
    atp_linked_ocr = {}
    for well in rate['Well'].unique():
        atp_linked_ocr[well] = ocr_pmol[well] - ocr_oligo_pmol[well]

    # Similarly, we can calculate proton leak OCR by subtracting OCR_rot from OCR_oligo
    proton_leak_ocr = {}
    for well in rate['Well'].unique():
        proton_leak_ocr[well] = ocr_oligo_pmol[well] - ocr_rot_pmol[well]

    # Finally, we can calculate ATP production rate from oxidative phosphorylation and glycolysis
    # ATP production from oxidative phosphorylation is calculated as:
    # ATP_oxphos = (ATP-linked OCR + (0.1 * OCR leak)) * P/O ratio * 2
    # Where P/O ratio is assumed to be 2.73 for complex I substrates
    p_o_ratio = 2.73
    atp_oxphos = {}
    for well in rate['Well'].unique():
        atp_oxphos[well] = (atp_linked_ocr[well] + (0.1 * proton_leak_ocr[well])) * p_o_ratio * 2  # in pmol ATP/min
    # ATP production from glycolysis is calculated as:
    # ATP_glycolysis = (Proton production rate from glycolysis - (0.38 * basal_OCR)) * 1.53
    atp_glycolysis = {}
    for well in rate['Well'].unique():
        atp_glycolysis[well] = (ecar_pmol[well] - (0.38 * basal_ocr[well])) * 1.53  # in pmol ATP/min
    # Total ATP production rate is the sum of ATP_oxphos and ATP_glycolysis
    total_atp_production = {}
    for well in rate['Well'].unique(): 
        total_atp_production[well] = atp_oxphos[well] + atp_glycolysis[well]  # in pmol ATP/min

    # Compile results into a DataFrame
    atp_results = pd.DataFrame({
        'Well': rate['Well'].unique(),
        "Group": [rate.loc[rate['Well'] == well, 'Group'].values[0] for well in rate['Well'].unique()],
        'Proton_Production_Rate_pmol_H_per_min': [ecar_pmol[well] for well in rate['Well'].unique()],
        'Basal_OCR_pmol_O2_per_min': [basal_ocr[well] for well in rate['Well'].unique()],
        'ATP_Linked_OCR_pmol_O2_per_min': [atp_linked_ocr[well] for well in rate['Well'].unique()],
        'Proton_Leak_OCR_pmol_O2_per_min': [proton_leak_ocr[well] for well in rate['Well'].unique()],
        'ATP_Production_OxPhos_pmol_ATP_per_min': [atp_oxphos[well] for well in rate['Well'].unique()],
        'ATP_Production_Glycolysis_pmol_ATP_per_min': [atp_glycolysis[well] for well in rate['Well'].unique()],
        'Total_ATP_Production_pmol_ATP_per_min': [total_atp_production[well] for well in rate['Well'].unique()]
    })
    atp_results.to_csv(os.path.join(res_dir, "ATP_production_results.csv"), index=False)
    print(f"\nATP production results saved to CSV ATP_production_results.csv in {res_dir}.")

    group_stats = atp_results.groupby('Group')[['ATP_Production_OxPhos_pmol_ATP_per_min', 'ATP_Production_Glycolysis_pmol_ATP_per_min']].mean()
    # Rename columns for clearer plot labels
    group_stats = group_stats.rename(columns={"ATP_Production_OxPhos_pmol_ATP_per_min": "OXPHOS",
                                                    "ATP_Production_Glycolysis_pmol_ATP_per_min": "Glycolysis"})

    # Plot stacked bar
    ax = group_stats.plot(kind='bar', stacked=True, figsize=(9,6), color=['#1f77b4', '#ff7f0e'])
    ax.set_ylabel('ATP production (pmol ATP/min)')
    ax.set_xlabel('Experimental Group')
    ax.set_title('ATP production per experimental group (mean across wells)')
    plt.tight_layout()
    # Save figure to results directory
    plt.savefig(os.path.join(res_dir, 'ATP_production_stacked_by_group.png'))
    print(f"ATP production stacked bar plot saved to ATP_production_stacked_by_group.png in {res_dir}.")
    label.configure(text="Processing complete. Results saved.")

def browsefiles():
    filename = filedialog.askopenfilename(initialdir = "/",
                                            title = "Select a File",
                                            filetypes = (("Excel files",
                                                            "*.xlsx"),
                                                         ("all files",
                                                            "*.*")))
    label_file_explorer.configure(text="File Opened: "+filename)
    print("File opened:", filename)
def getsavepath():
    save_path = filedialog.askdirectory()
    label_file_save.configure(text="Save Path: "+save_path)
    print("Save path selected:", save_path)

label_file_explorer = tk.Label(root,
                            text = "Select a File",
                            width = 25, height = 4,
                            fg = "blue")
button_explore = tk.Button(root,
                        text = "Browse Files",
                        command = browsefiles)
label_file_explorer.grid(column= 1, row = 1)
button_explore.grid(column = 1, row = 2)

label_file_save = tk.Label(root,
                            text = "Where to Save Results",
                            width = 25, height = 4,
                            fg = "blue")
button_save = tk.Button(root,
                        text = "Select Folder",
                        command = getsavepath)
label_file_save.grid(column= 1, row = 3)
button_save.grid(column = 1, row = 4)

oligo_entry = tk.Entry(root, width=5)
oligo_entry.grid(column=1, row=5)
oligo_label = tk.Label(root, text="Oligomycin injection number:")
oligo_label.grid(column=0, row=5)

rot_entry = tk.Entry(root, width=5)
rot_entry.grid(column=1, row=6)
rot_label = tk.Label(root, text="Rotenone/Antimycin A injection number:")
rot_label.grid(column=0, row=6)

button_run = tk.Button(root, text="Run", command=lambda: atp_prod_calc(label_file_explorer.cget("text").replace("File Opened: ",""), label_file_save.cget("text").replace("Save Path: ",""), int(oligo_entry.get()), int(rot_entry.get())))
button_run.grid(column= 1, row = 7)

button_exit = tk.Button(root, text="Exit", command=root.quit)
button_exit.grid(column= 1, row = 8)

root.mainloop()

