import matplotlib.pyplot as plt
import pandas as pd
import os,sys

dir_macs3=input('\nEnter directory where narrowpeak bed file output from macs3 is located.\n Do not put "\" at end:\n')
os.system('cp -R '+dir_macs3+'-filtered')
dir_files=dir_macs3+'_filtered'
file_list=os.listdir(dir_files)
print(file_list)

dir_plot=dir_files+'/plots'
os.system('mkdir '+dir_plot)
for file in file_list:
    if '.bed' in file:
         # Read BED file
        df = pd.read_csv(
            file,
            sep="\t",
            header=None
        )

        # Column 7 = signalValue
        signal = pd.to_numeric(
            df.iloc[:, 6],
            errors="coerce"
        )

        # Remove rows where column 7 is not numeric
        valid = signal.notna()
        df = df.loc[valid].copy()
        signal = signal.loc[valid]

        # -------------------------------------------------
        # Calculate 15th percentile cutoff
        # -------------------------------------------------

        cutoff = signal.quantile(0.15)

        keep = signal >= cutoff
        remove = signal < cutoff

        df_keep = df.loc[keep].copy()
        df_remove = df.loc[remove].copy()

        # -------------------------------------------------
        # Print statistics
        # -------------------------------------------------

        print("\n" + os.path.basename(file))
        print(f"15th percentile cutoff: {cutoff:.3f}")
        print(f"Total peaks: {len(signal)}")
        print(f"Removed: {remove.sum()} ({remove.mean()*100:.2f}%)")
        print(f"Kept: {keep.sum()} ({keep.mean()*100:.2f}%)")

        # -------------------------------------------------
        # Plot
        # -------------------------------------------------

        basename = os.path.basename(file)
        name = os.path.splitext(basename)[0]
        condition = basename.replace("peak.", "").replace(".narrowPeak.bed", "")

        
        fig, ax = plt.subplots(figsize=(8, 5))

        # All peaks
        ax.hist(
            signal,
            bins=30,
            alpha=0.6,
            label="All peaks"
        )

        # Removed peaks
        ax.hist(
            signal[remove],
            bins=30,
            alpha=0.8,
            label="Removed (<15th percentile)"
        )

        # Cutoff
        ax.axvline(
            cutoff,
            linestyle="--",
            linewidth=2,
            label=f"15th percentile = {cutoff:.2f}"
        )

        ax.set_title(name)
        ax.set_xlabel("MACS3 signalValue (column 7)")
        ax.set_ylabel("Number of peaks")
        ax.legend()

        plt.tight_layout()

        plot_file = os.path.join(
            dir_plot,
            f"peak.{condition}_filtered.png"
        )

        plt.savefig(
            plot_file,
            dpi=300
        )

        plt.close()

        # -------------------------------------------------
        # Save filtered BED
        # -------------------------------------------------

        filtered_file = os.path.join(
            dir_files,
            f"peak.{condition}_filtered.narrowPeak.bed"
        )

        df_keep.to_csv(
            filtered_file,
            sep="\t",
            header=False,
            index=False
        )
        os.system('rm '+file)
        print("Plot:", plot_file)
        print("Filtered BED files that can be used downstream:", filtered_file)

