import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import os
import argparse
import ast  # To safely parse the IonSeries string into a dictionary

def read_input_data(csv_file):
    """Read input CSV file into pandas DataFrame"""
    return pd.read_csv(csv_file)

def create_output_directory(output_dir):
    """Create output directory if it doesn't exist"""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

def define_colors():
    """Define color mapping for different ion types"""
    return {
        # Primary HCD ions (peptide backbone)
        "b": "#0066CC",  # Strong blue
        "y": "#00AA00",  # Strong green
        
        # Secondary HCD ions (diagnostic)
        "immonium": "#9370DB",  # Light purple
        "internal": "#DEB887",  # Light brown
        "oxonium": "#FF69B4",  # Light pink
        
        # Glycan-related ions
        "Y": "#FF0000",  # Strong red
        "Y_peptide_loss": "#FFB6C1",  # Light red
        "B": "#FFA500",  # Orange
        
        # ETD ions
        "c": "#FFD700",  # Gold
        "z": "#FF1493",  # Deep pink
    }

def define_ion_properties():
    """Define properties for each ion type"""
    return {
        # HCD ions configuration
        'hcd': {
            # Primary ions (full opacity)
            'b': {'length': 40, 'style': '-', 'alpha': 1.0, 'zorder': 3},
            'y': {'length': 35, 'style': '-', 'alpha': 1.0, 'zorder': 3},
            
            # Secondary ions (reduced opacity)
            'immonium': {'length': 15, 'style': ':', 'alpha': 0.6, 'zorder': 1},
            'internal': {'length': 20, 'style': ':', 'alpha': 0.6, 'zorder': 1},
            'oxonium': {'length': 45, 'style': ':', 'alpha': 0.7, 'zorder': 2},
            
            # Glycan-related ions
            'Y': {'length': 30, 'style': '--', 'alpha': 0.9, 'zorder': 2},
            'Y_peptide_loss': {'length': 25, 'style': ':', 'alpha': 0.5, 'zorder': 1},
            'B': {'length': 40, 'style': '-.', 'alpha': 0.8, 'zorder': 2},
        },
        
        # ETD ions configuration
        'etd': {
            'c': {'length': 40, 'style': '-', 'alpha': 1.0, 'zorder': 3},
            'z': {'length': 35, 'style': '-', 'alpha': 1.0, 'zorder': 3},
            'Y': {'length': 30, 'style': '--', 'alpha': 0.9, 'zorder': 2},
            'Y_peptide_loss': {'length': 25, 'style': ':', 'alpha': 0.5, 'zorder': 1},
        }
    }

def create_mass_spectrum_plots(df_row, color_mapping, base_title):
    """Create separate HCD and ETD mass spectrum plots"""
    
    # Create two figures with white background
    fig_hcd, ax_hcd = plt.subplots(figsize=(15, 8), facecolor='white')
    fig_etd, ax_etd = plt.subplots(figsize=(15, 8), facecolor='white')
    ax_hcd.set_facecolor('white')
    ax_etd.set_facecolor('white')
    
    ion_series = df_row['IonSeries'].iloc[0]
    ion_properties = define_ion_properties()
    
    # Plot HCD spectrum
    legend_handles_hcd = {}
    for ion_type, values in ion_series.items():
        if ion_type in ion_properties['hcd']:
            _plot_ion_series(ax_hcd, values, ion_type, ion_properties['hcd'][ion_type], 
                           color_mapping, legend_handles_hcd)
    
    # Plot ETD spectrum
    legend_handles_etd = {}
    for ion_type, values in ion_series.items():
        if ion_type in ion_properties['etd']:
            _plot_ion_series(ax_etd, values, ion_type, ion_properties['etd'][ion_type], 
                           color_mapping, legend_handles_etd)
    
    # Configure HCD plot
    _configure_plot(ax_hcd, legend_handles_hcd, base_title + "\nHCD Fragmentation", "HCD Ion Types")
    
    # Configure ETD plot
    _configure_plot(ax_etd, legend_handles_etd, base_title + "\nETD Fragmentation", "ETD Ion Types")
    
    plt.tight_layout()
    return fig_hcd, fig_etd

def _plot_ion_series(ax, values, ion_type, props, color_mapping, legend_handles):
    """Plot a series of ions with consistent styling"""
    color = color_mapping.get(ion_type, "black")
    
    if isinstance(values, list):
        for i, mz in enumerate(values):
            ion_number = i + 1
            ion_label = f"{ion_type}{ion_number}"
            _plot_single_ion(ax, mz, ion_label, ion_type, props, color, legend_handles)
    elif isinstance(values, dict):
        for ion_name, mz in values.items():
            _plot_single_ion(ax, mz, ion_name, ion_type, props, color, legend_handles)

def _plot_single_ion(ax, mz, ion_label, ion_type, props, color, legend_handles):
    """Plot a single ion with specified properties"""
    # Plot vertical line with specified properties
    ax.vlines(mz, 0, props['length'], color=color, lw=2, 
             linestyle=props['style'], alpha=props['alpha'], 
             zorder=props['zorder'])
    
    # Add text label with white outline
    text = ax.text(mz, props['length'], ion_label, rotation=90,
                  verticalalignment='bottom', horizontalalignment='center',
                  fontsize=8, color=color, alpha=props['alpha'],
                  zorder=props['zorder'])
    
    # Add white outline to text for better visibility
    text.set_path_effects([
        path_effects.Stroke(linewidth=3, foreground='white'),
        path_effects.Normal()
    ])
    
    # Add to legend if not already present
    if ion_type not in legend_handles:
        legend_handles[ion_type] = plt.Line2D(
            [0], [0], color=color, lw=2, 
            linestyle=props['style'], 
            alpha=props['alpha'],
            label=_get_legend_label(ion_type)
        )

def _get_legend_label(ion_type):
    """Generate descriptive legend labels"""
    labels = {
        'b': 'b-ions (HCD backbone)',
        'y': 'y-ions (HCD backbone)',
        'c': 'c-ions (ETD backbone)',
        'z': 'z-ions (ETD backbone)',
        'Y': 'Y-ions (glycopeptide)',
        'Y_peptide_loss': 'Y-ions with peptide loss',
        'B': 'B-ions (glycan)',
        'immonium': 'Immonium ions',
        'internal': 'Internal fragments',
        'oxonium': 'Oxonium ions'
    }
    return labels.get(ion_type, ion_type)

def _configure_plot(ax, legend_handles, title, legend_title):
    """Configure plot appearance"""
    # Sort legend handles by ion type importance
    ion_order = ['b', 'y', 'c', 'z', 'Y', 'B', 'oxonium', 'immonium', 'internal', 'Y_peptide_loss']
    sorted_handles = sorted(
        legend_handles.items(),
        key=lambda x: ion_order.index(x[0]) if x[0] in ion_order else len(ion_order)
    )
    
    # Add legend with custom ordering
    ax.legend(
        [h for _, h in sorted_handles],
        [h.get_label() for _, h in sorted_handles],
        title=legend_title,
        loc="upper right",
        fontsize=8,
        title_fontsize=10,
        bbox_to_anchor=(1.15, 1),
        borderaxespad=0.
    )
    
    # Configure axes
    ax.set_xlabel("Calculated m/z values (Da)", fontsize=10)
    ax.set_ylabel("Relative Intensity (%)", fontsize=10)
    ax.set_title(title, fontsize=12, pad=20)
    ax.set_ylim(0, 100)
    
    # Add grid for better readability
    ax.grid(True, linestyle=':', alpha=0.3, zorder=0)
    
    # Add spines
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.5)

def save_plot(fig, output_file):
    """Save plot to output location with high DPI"""
    fig.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"Plot saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Plot separate HCD and ETD mock mass spectra from glycopeptide ion series CSV file.")
    parser.add_argument('-i', '--input', required=True, help="Input CSV file containing glycopeptide ion series info.")
    parser.add_argument('-o', '--output', default="mock_mass_spectra", help="Output directory to save the plots.")
    args = parser.parse_args()
    
    # Create output directory
    create_output_directory(args.output)
    
    # Read input data
    df = read_input_data(args.input)
    
    # Check required columns
    required_columns = ['IonSeries', 'ProteinID', 'Peptide', 'Composition', 'GlyToucan_AC',
                       'GlycopeptideMass', 'GlycanMass', 'PeptideMass', 'Site']
    for col in required_columns:
        if col not in df.columns:
            raise KeyError(f"Missing required column: {col}")
    
    # Parse IonSeries column
    df['IonSeries'] = df['IonSeries'].apply(ast.literal_eval)
    
    # Define color mapping
    color_mapping = define_colors()
    
    # Process each row
    for index, row in df.iterrows():
        # Skip long peptides
        if len(row['Peptide']) > 50:
            continue
        
        # Create base title
        base_title = (
            f"Glycopeptide Sequence Finder: Mock Mass Spectrum\n"
            f"Protein: {row['ProteinID']}, Site: {int(row['Site'])}\n"
            f"Peptide: {row['Peptide']}, Mass: {row['PeptideMass']:.2f} Da\n"
            f"Glycan: {row['Composition']}, Mass: {row['GlycanMass']:.2f} Da"
        )
        
        # Create single row dataframe for plotting
        row_df = pd.DataFrame([row])
        
        # Generate plots
        fig_hcd, fig_etd = create_mass_spectrum_plots(row_df, color_mapping, base_title)
        
        # Create output filenames
        base_filename = f"{row['ProteinID']}_{int(row['Site'])}_{row['Peptide']}_{row['GlyToucan_AC']}".replace('|', '_')
        hcd_output = os.path.join(args.output, f"{base_filename}_HCD_spectrum.png")
        etd_output = os.path.join(args.output, f"{base_filename}_ETD_spectrum.png")
        
        # Save plots
        save_plot(fig_hcd, hcd_output)
        save_plot(fig_etd, etd_output)

if __name__ == "__main__":
    main()