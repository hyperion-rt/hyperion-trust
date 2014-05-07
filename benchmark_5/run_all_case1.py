from hyperion.model import Model
from cli            import command_line_arguments, handle_command_line_arguments, filename
from plot_results   import plot_results
from setup_model    import setup_model
from export_to_fits import export_to_fits

cli = command_line_arguments()

# Print comments to screen:
cli.verbose             = 0

# Resolution in pc:
cli.resolution          = 100.0

# Photons per grid cell:
cli.photons_temperature = 3
cli.photons_raytracing  = 3
cli.photons_imaging     = 3

cli.case                = 1

for cli.opticaldepth in [0.256434, 1.28217, 6.41084]:
    for cli.mode in ["temperature", "images", "seds"]:
        
        # Setup step:
        setup_model(cli);
        
        # Compute step:
        file = filename(cli, cli.mode)
        model = Model.read(file+".rtin")
        model.write("temp.rtin")
        model.run(filename=file+".rtout", logfile=file+".hyout", overwrite=True)
        
        # Plotting step:
        #TODO: if(cli.mode == "temperature"): visualize_setup + visualize_results => plot_setup
        if(cli.mode != "temperature"):
            plot_results(cli);
        
        # Export to fits format:
        if(cli.mode == "images"):
            export_to_fits(cli);
