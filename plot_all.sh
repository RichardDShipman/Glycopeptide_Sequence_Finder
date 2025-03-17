
input_dir="digested_glycopeptide_library"
cores=4

ls ${input_dir}/*.csv | xargs -I {} -P ${cores} python plot_mock_mass_spectra.py -i "{}"

#